"""The optional AI fallback for `ask`.  No network: the transport is faked."""
import io
import json
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from prospector.core import ai, cache, http          # noqa: E402


def answer(steps=(), cannot="", stop="end_turn"):
    text = json.dumps({"steps": [{"command": c, "description": d} for c, d in steps],
                       "cannot": cannot})
    return {"model": "claude-opus-5-5", "stop_reason": stop,
            "content": [{"type": "thinking", "thinking": ""}, {"type": "text", "text": text}],
            "usage": {"input_tokens": 2100, "output_tokens": 60}}


class FakeService(object):
    def __init__(self, reply=None, error=None):
        self.reply, self.error, self.calls = reply, error, []

    def __call__(self, url, body, headers, timeout):
        self.calls.append((url, json.loads(body.decode("utf-8")), headers, timeout))
        if self.error:
            raise self.error
        return 200, json.dumps(self.reply).encode("utf-8"), {}


class Base(unittest.TestCase):
    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        self.saved_env = os.environ.pop("ANTHROPIC_API_KEY", None)
        self.saved_transport = http.transport

    def tearDown(self):
        http.transport = self.saved_transport
        if self.saved_env is not None:
            os.environ["ANTHROPIC_API_KEY"] = self.saved_env
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)


class CoreTests(Base):
    def test_off_until_a_key_is_stored(self):
        self.assertFalse(ai.enabled())
        with self.assertRaises(ai.AIError) as caught:
            ai.ask("hello", "ref", "ctx")
        self.assertEqual(caught.exception.kind, "nokey")
        ai.set_key("  sk-ant-test-1234567890  \n")
        self.assertTrue(ai.enabled())
        self.assertEqual(ai.api_key(), "sk-ant-test-1234567890")
        self.assertEqual(ai.masked(ai.api_key()), "...7890")
        if os.name == "posix":
            self.assertEqual(os.stat(ai._key_path()).st_mode & 0o777, 0o600)
        ai.set_enabled(False)
        self.assertFalse(ai.enabled())
        ai.set_enabled(True)
        self.assertTrue(ai.forget_key())
        self.assertFalse(ai.enabled())

    def test_the_request_follows_the_documented_shape(self):
        ai.set_key("sk-ant-test-1234567890")
        http.transport = service = FakeService(answer([("viewas sticks", "draw as sticks")]))
        plan = ai.ask("mostrami la proteina a bastoncini", "REFERENCE", "object 1UBQ")
        self.assertEqual(plan.steps, [("viewas sticks", "draw as sticks")])
        url, body, headers, timeout = service.calls[0]
        self.assertEqual(url, "https://api.anthropic.com/v1/messages")
        self.assertEqual(headers["x-api-key"], "sk-ant-test-1234567890")
        self.assertEqual(headers["anthropic-version"], "2026-10-21")
        self.assertEqual(body["model"], "claude-opus-5-5")
        self.assertEqual(body["output_config"]["format"]["type"], "json_schema")
        self.assertEqual(body["output_config"]["effort"], "low")
        self.assertEqual(body["fallbacks"], "default")
        self.assertEqual(headers["anthropic-beta"], "server-side-fallback-2026-07-01")
        for removed in ("temperature", "top_p", "top_k", "thinking", "tool_choice"):
            self.assertNotIn(removed, body)
        self.assertIn("REFERENCE", body["system"][0]["text"])
        self.assertEqual(body["system"][0]["cache_control"], {"type": "ephemeral"})
        # what varies is in the user message, so the system prompt stays cacheable
        self.assertIn("object 1UBQ", body["messages"][0]["content"])
        self.assertIn("bastoncini", body["messages"][0]["content"])
        self.assertNotIn("1UBQ", body["system"][0]["text"])
        self.assertEqual(len(service.calls), 1)             # one attempt, no retry

    def test_other_models_get_no_fallback_parameter(self):
        ai.set_key("sk-ant-test-1234567890")
        ai.set_model("claude-haiku-5-5")
        http.transport = service = FakeService(answer([("status", "x")]))
        ai.ask("hm", "ref", "ctx")
        _, body, headers, _ = service.calls[0]
        self.assertEqual(body["model"], "claude-haiku-5-5")
        self.assertNotIn("fallbacks", body)
        self.assertNotIn("anthropic-beta", headers)
        with self.assertRaises(ai.AIError):
            ai.set_model("gpt-something; rm -rf")

    def test_problems_become_messages_a_student_can_act_on(self):
        ai.set_key("sk-ant-test-1234567890")
        cases = [(http.NetworkError("x", status=401), "key", "prospector ai key"),
                 (http.NetworkError("x", status=429), "busy", "minute"),
                 (http.NetworkError("x", status=529), "busy", "overloaded"),
                 (http.NetworkError("x", offline=True), "offline", "internet"),
                 (http.NetworkError("did not answer"), "offline", "in time"),
                 (http.NetworkError("credit balance is too low", status=400), "credit",
                  "credit")]
        for error, kind, words in cases:
            http.transport = FakeService(error=error)
            with self.assertRaises(ai.AIError) as caught:
                ai.ask("hm", "ref", "ctx")
            self.assertEqual(caught.exception.kind, kind)
            self.assertIn(words, str(caught.exception))
        http.transport = FakeService(answer(stop="refusal"))
        with self.assertRaises(ai.AIError) as caught:
            ai.ask("hm", "ref", "ctx")
        self.assertEqual(caught.exception.kind, "refused")
        http.transport = FakeService({"stop_reason": "end_turn",
                                      "content": [{"type": "text", "text": "not json"}]})
        with self.assertRaises(ai.AIError):
            ai.ask("hm", "ref", "ctx")

    def test_an_unsupported_extra_is_dropped_and_the_question_resent(self):
        ai.set_key("sk-ant-test-1234567890")
        calls = []

        def picky(url, body, headers, timeout):
            sent = json.loads(body.decode("utf-8"))
            calls.append(sent)
            if "fallbacks" in sent:
                raise http.NetworkError("HTTP 400. fallbacks: unknown beta", status=400)
            if "format" in sent["output_config"]:
                raise http.NetworkError("HTTP 400. output_config.format: json_schema is "
                                        "not supported", status=400)
            text = "```json\n%s\n```" % json.dumps(
                {"steps": [{"command": "status", "description": "x"}], "cannot": ""})
            return 200, json.dumps({"stop_reason": "end_turn", "content": [
                {"type": "text", "text": text}]}).encode("utf-8"), {}

        http.transport = picky
        plan = ai.ask("hm", "ref", "ctx")
        self.assertEqual(plan.steps, [("status", "x")])
        self.assertEqual(len(calls), 3)
        self.assertIn("Reply with one JSON object", calls[2]["system"][0]["text"])

    def test_only_safe_commands_pass_the_check(self):
        allowed = ["focus active site of p53", "whatif p53 R248Q", "viewas sticks",
                   "find lysozyme res<2.0 organism=human", "color red, ss h",
                   "hide everything, solvent", "set cartoon_transparency, 0.5",
                   "distance hb, polymer, organic, mode=2", "zoom chain A", "clear",
                   "session save my_work", "show sticks, resn HIS and chain A"]
        for command in allowed:
            self.assertIsNone(ai.check(command), command)
        refused = ["delete all", "run /tmp/evil.py", "python", "save /etc/passwd",
                   "viewas sticks; delete all", "focus all\ndelete all", "fetch 1ubq",
                   "label all, __import__('os').system('x')", "alter all, b=0",
                   "set fetch_path, /tmp", "session delete my_work", "fav rm 1ubq",
                   "cmd.delete('all')", "@script.pml", "batch list.txt", "ask clear",
                   "prospector ai off", "reinitialize", "", "x" * 300, "system ls",
                   "color red, `id`", "load http://example.org/x.pse", "cd /", "quit"]
        for command in refused:
            self.assertIsNotNone(ai.check(command), command)
        plan = ai.Plan([("viewas sticks", "a"), ("delete all", "b")])
        good, bad = ai.validate(plan)
        self.assertEqual([c for c, _ in good], ["viewas sticks"])
        self.assertEqual(bad[0][0], "delete all")


class AskTests(Base):
    def setUp(self):
        Base.setUp(self)
        from pymol import cmd
        from prospector.pymol_layer import aiask, ask, commands
        self.cmd, self.ask, self.commands, self.aiask = cmd, ask, commands, aiask
        commands.register()
        cmd.delete("all")
        cmd.fab("MKHACDEFGHIKLMNPQRSTVWYHCAAC", "pep", ss=1)
        commands._remember("mypep", "pep")

    def tearDown(self):
        self.cmd.delete("all")
        Base.tearDown(self)

    def run_line(self, line):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            self.cmd.do(line, echo=0)
        return buffer.getvalue()

    def test_without_a_key_nothing_is_sent(self):
        http.transport = service = FakeService(answer([("viewas sticks", "x")]))
        out = self.run_line("ask florp the glorbles")
        self.assertEqual(service.calls, [])
        self.assertIn("I'm not sure", out)
        self.assertIn("prospector ai", out)

    def test_the_rules_go_first_and_the_ai_is_not_asked(self):
        ai.set_key("sk-ant-test-1234567890")
        http.transport = service = FakeService(answer([("clear", "x")]))
        out = self.run_line("ask show it as sticks")
        self.assertEqual(service.calls, [])
        self.assertIn("viewas sticks", out)
        self.assertNotIn("(AI)", out)

    def test_the_ai_is_asked_when_the_rules_give_up(self):
        ai.set_key("sk-ant-test-1234567890")
        http.transport = service = FakeService(answer(
            [("viewas mypep spheres", "draw mypep as spheres"),
             ("color blue, ss h", "colour the helices blue")]))
        out = self.run_line("ask florp the glorbles of mypep")
        self.assertEqual(len(service.calls), 1)
        self.assertIn("Asking the AI", out)
        self.assertIn("I UNDERSTOOD (AI)  draw mypep as spheres", out)
        self.assertIn("RUNNING     viewas mypep spheres", out)
        self.assertIn("RUNNING     color blue, ss h", out)
        self.assertTrue(self.cmd.count_atoms("pep and rep spheres") > 0)
        _, body, _, _ = service.calls[0]
        message = body["messages"][0]["content"]
        self.assertIn("florp the glorbles of mypep", message)
        self.assertIn("object pep", message)            # what is loaded
        self.assertIn("mypep", message)                 # the name the user gave it
        self.assertIn("whatif", body["system"][0]["text"])   # the command manual
        self.assertIn("scout", body["system"][0]["text"])

    def test_a_forbidden_answer_runs_nothing(self):
        ai.set_key("sk-ant-test-1234567890")
        http.transport = FakeService(answer([("viewas spheres", "ok"),
                                             ("delete all", "remove everything")]))
        out = self.run_line("ask florp the glorbles")
        self.assertIn("AI NOT USED", out)
        self.assertIn("delete all", out)
        self.assertIn("pep", self.cmd.get_names("objects"))
        self.assertEqual(self.cmd.count_atoms("pep and rep spheres"), 0)   # all or nothing

    def test_cannot_and_errors_are_reported_plainly(self):
        ai.set_key("sk-ant-test-1234567890")
        http.transport = FakeService(answer(cannot="Prospector cannot simulate folding."))
        out = self.run_line("ask florp the glorbles")
        self.assertIn("THE AI SAYS  Prospector cannot simulate folding.", out)
        self.assertIn("Nothing was run", out)
        http.transport = FakeService(error=http.NetworkError("x", offline=True))
        out = self.run_line("ask florp the glorbles")
        self.assertIn("AI NOT USED", out)
        self.assertIn("THINGS I UNDERSTAND", out)

    def test_how_do_i_only_shows_the_ai_commands(self):
        ai.set_key("sk-ant-test-1234567890")
        http.transport = FakeService(answer([("viewas spheres", "draw as spheres")]))
        out = self.run_line("ask how do I florp the glorbles")
        self.assertIn("TYPE        viewas spheres", out)
        self.assertIn("Nothing was run", out)
        self.assertEqual(self.cmd.count_atoms("pep and rep spheres"), 0)

    def test_askai_skips_the_rules(self):
        # without a key it explains how to get one and sends nothing
        http.transport = service = FakeService(answer([("viewas spheres", "x")]))
        out = self.run_line("askai show it as sticks")
        self.assertEqual(service.calls, [])
        self.assertIn("needs an AI key", out)
        self.assertIn("ask show it as sticks", out)
        # with a key even a sentence the rules know goes to the AI
        ai.set_key("sk-ant-test-1234567890")
        out = self.run_line("askai show it as sticks")
        self.assertEqual(len(service.calls), 1)
        self.assertIn("I UNDERSTOOD (AI)", out)
        self.assertIn("RUNNING     viewas spheres", out)     # the AI's answer, not the rule's
        self.assertNotIn("built-in rules did not understand", out)
        self.assertIn("show it as sticks", service.calls[0][1]["messages"][0]["content"])
        # errors are reported and point back to ask
        http.transport = FakeService(error=http.NetworkError("x", status=429))
        out = self.run_line("askai florp")
        self.assertIn("AI NOT USED", out)
        self.assertIn("MEANWHILE    ask florp", out)
        self.assertIn("straight to the AI", self.run_line("askai"))
        self.assertIn("askai", self.run_line("info"))

    def test_an_ai_answer_is_remembered_and_reused_without_the_ai(self):
        ai.set_key("sk-ant-test-1234567890")
        http.transport = service = FakeService(answer(
            [("viewas mypep spheres", "draw mypep as spheres")]))
        first = self.run_line("ask florp the glorbles of mypep")
        self.assertEqual(len(service.calls), 1)
        self.assertIn("Remembered", first)
        second = self.run_line("ask Florp the glorbles of mypep!")
        self.assertEqual(len(service.calls), 1)              # no second request
        self.assertIn("I UNDERSTOOD (learned)  draw mypep as spheres", second)
        self.assertIn("RUNNING     viewas mypep spheres", second)
        # it works with the AI switched off, and with no key at all
        ai.forget_key()
        self.assertIn("(learned)", self.run_line("ask florp the glorbles of mypep"))
        listing = self.run_line("prospector ai learned")
        self.assertIn("florp the glorbles of mypep", listing)
        self.assertIn("Forgot 1 learned sentence", self.run_line("prospector ai unlearn"))
        self.assertIn("I'm not sure", self.run_line("ask florp the glorbles of mypep"))

    def test_a_learned_answer_follows_the_structure_being_worked_on(self):
        """The real case: the AI names the loaded structure although the
        sentence does not.  The answer is stored for "the current structure"."""
        ai.set_key("sk-ant-test-1234567890")
        http.transport = service = FakeService(answer(
            [("color grey, pep", "Colora di grigio tutta la struttura."),
             ("color blue, pep and ss h", "Colora di blu le eliche alfa.")]))
        sentence = "ask colora di blu le eliche e di grigio tutto il resto"
        self.assertIn("Remembered", self.run_line(sentence))
        stored = list(ai.learned_all().values())[0]["steps"]
        self.assertEqual(stored[1][0], "color blue, {current} and ss h")
        again = self.run_line(sentence)
        self.assertEqual(len(service.calls), 1)
        self.assertIn("(learned)", again)
        self.assertIn("RUNNING     color blue, pep and ss h", again)
        # another structure is loaded later: the answer follows it
        self.cmd.fab("ACDEFGHIKLMNPQRSTVWY", "other", ss=1)
        self.commands._last[0] = "other"
        self.assertIn("RUNNING     color blue, other and ss h", self.run_line(sentence))
        # nothing loaded: the learned answer cannot apply, so the AI is asked
        self.cmd.delete("all")
        self.commands._last[0] = None
        self.run_line(sentence)
        self.assertEqual(len(service.calls), 2)
        # a sentence that names the structure keeps the name
        self.assertEqual(ai.generalise("colour pep red", [("color red, pep", "x")], "pep"),
                         [("color red, pep", "x")])

    def test_answers_that_depend_on_the_screen_are_not_remembered(self):
        # points back to an earlier step
        self.assertFalse(ai.learnable("the same for chain B", [("color red, chain B", "x")]))
        # names a loaded object the sentence never mentions
        self.assertFalse(ai.learnable("colour the helices", [("color red, pep and ss h", "x")],
                                      ["pep"]))
        self.assertTrue(ai.learnable("colour the helices of pep red",
                                     [("color red, pep and ss h", "x")], ["pep"]))
        self.assertTrue(ai.learnable("colour the helices red", [("color red, ss h", "x")],
                                     ["pep"]))
        # a stored command is checked again before it is reused
        ai.learned_put("do the thing", [("delete all", "x")])
        self.assertIsNone(ai.learned_get("do the thing"))

    def test_switched_off_means_nothing_is_sent(self):
        ai.set_key("sk-ant-test-1234567890")
        out = self.run_line("prospector ai off")
        self.assertIn("OFF", out)
        http.transport = service = FakeService(answer([("clear", "x")]))
        self.run_line("ask florp the glorbles")
        self.assertEqual(service.calls, [])
        self.assertIn("ON", self.run_line("prospector ai on"))

    def test_settings_never_show_or_accept_the_key_on_the_command_line(self):
        ai.set_key("sk-ant-secret-ABCDEFGH1234")
        out = self.run_line("prospector ai")
        self.assertIn("...1234", out)
        self.assertNotIn("sk-ant-secret", out)
        self.assertIn("WHAT IS SENT", out)
        out = self.run_line("prospector ai key sk-ant-typed-in-the-log")
        self.assertIn("Do not type the key", out)
        self.assertEqual(ai.api_key(), "sk-ant-secret-ABCDEFGH1234")
        http.transport = FakeService(answer([("viewas sticks", "x")]))
        self.assertIn("WORKING", self.run_line("prospector ai test"))
        http.transport = FakeService(error=http.NetworkError("x", status=401))
        self.assertIn("NOT WORKING", self.run_line("prospector ai test"))

    def test_the_manual_is_stable_and_complete(self):
        first = self.aiask.reference()
        self.aiask._reference[0] = None
        self.assertEqual(first, self.aiask.reference())      # same bytes: cacheable
        for name in ai.PROSPECTOR_COMMANDS:
            self.assertIn("\n%s - " % name, "\n" + first)
        self.assertLess(len(first), 40000)


if __name__ == "__main__":
    unittest.main()
