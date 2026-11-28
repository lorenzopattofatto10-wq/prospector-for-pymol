"""Run slow jobs (network) off PyMOL's GUI thread, and bring the results back.

* In the PyMOL GUI: jobs run one after another in a single background
  thread; each result is handed back to the GUI thread through a Qt signal,
  so PyMOL never freezes and output appears in the order commands were typed.
* Without a GUI event loop (``pymol -c``, scripts run from a thread, unit
  tests) jobs simply run in place.
"""
from __future__ import annotations

import queue
import threading
import traceback

_async_enabled = True
class _Line(object):
    """The waiting line of jobs.  Normally first come, first served; a command
    that expands into other commands (``ask``) can put those at the FRONT, so
    they run right after it and before anything typed later."""

    def __init__(self):
        self._items = []
        self._ready = threading.Condition()
        self.front = None        # position for front insertion, or None

    def put(self, item):
        with self._ready:
            if self.front is None:
                self._items.append(item)
            else:
                self._items.insert(self.front, item)
                self.front += 1
            self._ready.notify()

    def get(self):
        with self._ready:
            while not self._items:
                self._ready.wait()
            return self._items.pop(0)


_jobs = _Line()


class next_in_line(object):
    """``with next_in_line():`` - commands started inside run before the ones
    already waiting, in the order they were started."""

    def __enter__(self):
        _jobs.front = 0

    def __exit__(self, *exc):
        _jobs.front = None
        return False
_thread = None
_bridge = None
_pending = 0
_lock = threading.Lock()


def set_async(enabled):
    global _async_enabled
    _async_enabled = bool(enabled)


def is_async():
    return _async_enabled


def _qt():
    try:
        from pymol.Qt import QtCore, QtWidgets
        return QtCore, QtWidgets
    except Exception:
        return None, None


def _on_gui_thread():
    """True when a Qt event loop exists and we are being called from its thread."""
    QtCore, QtWidgets = _qt()
    if QtCore is None:
        return False
    app = QtWidgets.QApplication.instance()
    if app is None:
        return False
    return QtCore.QThread.currentThread() == app.thread()


def init():
    """Create the object that carries results back to the GUI thread.
    Call once from the GUI thread (plugin start-up)."""
    global _bridge
    QtCore, QtWidgets = _qt()
    if _bridge is not None or QtCore is None or QtWidgets.QApplication.instance() is None:
        return

    class Bridge(QtCore.QObject):
        deliver = QtCore.Signal(object)

        def __init__(self):
            QtCore.QObject.__init__(self)
            self.deliver.connect(self._run, QtCore.Qt.QueuedConnection)

        def _run(self, function):
            function()

    bridge = Bridge()
    app = QtWidgets.QApplication.instance()
    if bridge.thread() != app.thread():
        bridge.moveToThread(app.thread())
    _bridge = bridge


def _finish(done, result, error, finished=None):
    global _pending
    try:
        done(result, error)
    except Exception:
        print("Prospector: unexpected error while showing a result:\n" + traceback.format_exc())
    finally:
        with _lock:
            _pending -= 1
        if finished is not None:
            finished.set()      # the worker may now start the next command


def _deliver(done, result, error):
    """Hand a result to the GUI thread and WAIT until it has been shown.  One
    command finishes completely - including any follow-up commands it starts -
    before the next one begins, so output never interleaves or overtakes."""
    finished = threading.Event()
    _bridge.deliver.emit(lambda: _finish(done, result, error, finished))
    finished.wait()


def _still_working():
    print("Prospector is working on it (asking the online databases) ...")


def _run_before(before, box, ready):
    try:
        box["go"] = before() is not False
    except Exception:
        box["error"] = traceback.format_exc()
    finally:
        ready.set()


def _loop():
    while True:
        before, job, done = _jobs.get()
        result, error = None, None
        if before is not None:
            # `before` must see PyMOL as the previous command left it, so it runs
            # on the GUI thread, but only now that it is this command's turn.
            box, ready = {}, threading.Event()
            _bridge.deliver.emit(lambda b=before, x=box, r=ready: _run_before(b, x, r))
            ready.wait()
            if box.get("error"):
                _deliver(done, None, box["error"])
                continue
            if not box.get("go"):
                _deliver(lambda r, e: None, None, None)
                continue
        # if the job is slow (usually: waiting for an online database), say so
        # once, so the user knows something is happening
        waiting = threading.Timer(1.5, lambda: _bridge.deliver.emit(_still_working))
        waiting.daemon = True
        waiting.start()
        try:
            result = job()
        except Exception:
            error = traceback.format_exc()
        waiting.cancel()
        _deliver(done, result, error)


def submit(job, done, before=None):
    """Run ``job()``; later call ``done(result, error_text)`` on the calling side.
    ``error_text`` is None unless the job raised an unexpected exception.

    ``before`` (optional) runs on the GUI thread just before the job starts -
    after every earlier command has completely finished.  Commands that look
    at loaded objects use it, so `grab X` followed by `features X` works even
    when both lines are pasted at once.  If it returns False the job is skipped."""
    global _thread, _pending
    with _lock:
        _pending += 1
    if _async_enabled and _on_gui_thread():
        init()
        if _bridge is not None:
            if _thread is None or not _thread.is_alive():
                _thread = threading.Thread(target=_loop, name="prospector-worker")
                _thread.daemon = True
                _thread.start()
            _jobs.put((before, job, done))
            return True
    result, error = None, None
    try:
        if before is not None and before() is False:
            done = lambda r, e: None
        else:
            result = job()
    except Exception:
        error = traceback.format_exc()
    _finish(done, result, error)
    return False


def later(function, delay_ms=300):
    """Run ``function`` on the GUI thread a moment from now, after PyMOL's
    window has redrawn; without a window, at once."""
    if _bridge is None:
        function()
        return
    QtCore, _ = _qt()

    def schedule():
        QtCore.QTimer.singleShot(delay_ms, function)
    _bridge.deliver.emit(schedule)


def pending():
    return _pending
