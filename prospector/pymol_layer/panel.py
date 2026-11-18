"""A small Qt panel.  It only builds command lines and runs the same PyMOL
commands you could type yourself - and shows you the command, so you learn it."""
from __future__ import annotations

from pymol import cmd

from . import commands as base_commands
from pymol.Qt import QtCore, QtGui, QtWidgets

from ..core import rank

COLUMNS = ["PDB code", "Score /100", "Method", "Detail (A)", "Year", "Organism", "Chains",
           "Bound molecules", "Title"]
_CARDS = ("about", "find", "info", "sequence", "status", "start", "history", "fav",
          "session", "similar", "glossary", "prospector", "compare")


def _readable(command, lines, hidden):
    """The text of a finished command as a page that is easy to read: each
    labelled block becomes a short paragraph with its label in bold.  Cards
    and tables keep their fixed-width layout."""
    import html
    import re
    head = ("<p style='color:gray;margin:0'>You ran</p><p style='margin:0 0 8px 0'><tt><b>%s"
            "</b></tt></p>" % html.escape(command))
    more = ("<p style='color:gray'>(%d more lines - press 'Show all of it')</p>" % hidden
            if hidden else "")
    first = command.split()[0] if command.split() else ""
    if first in _CARDS:
        return head + "<pre>%s</pre>" % html.escape("\n".join(lines)) + more
    blocks = []
    for line in lines:
        if not line.strip() or set(line.strip()) <= set("-="):
            continue
        labelled = re.match(r"^([A-Z][A-Za-z ?]{1,15}?)(?::|\s{2,})\s*(\S.*)$", line)
        row = re.match(r"^\s{2,}([A-Z][A-Za-z\-]{1,11})\s{2,}(\S.*)$", line)
        if labelled:
            blocks.append([labelled.group(1).strip(), labelled.group(2).strip(), []])
        elif row and blocks:                     # "   Charge      positive -> none"
            blocks[-1][2].append((row.group(1), row.group(2).strip()))
        elif line.startswith(" ") and blocks:      # one more item of the same block
            if blocks[-1][2]:
                last = blocks[-1][2][-1]
                blocks[-1][2][-1] = (last[0], last[1] + " " + line.strip())
            else:
                blocks[-1][1] += "\n" + line.strip()
        else:
            blocks.append(["", line.strip(), []])
    parts = []
    for label, text, rows in blocks:
        text = html.escape(text).replace("\n", "<br>")
        if label:
            parts.append("<p style='margin:6px 0 0 0'><b>%s</b><br>%s</p>"
                         % (html.escape(label.capitalize()), text))
        else:
            parts.append("<p style='margin:6px 0 0 0'>%s</p>" % text)
        if rows:
            parts.append("<table cellspacing='0' cellpadding='2'>%s</table>" % "".join(
                "<tr><td style='color:gray'>%s&nbsp;&nbsp;</td><td>%s</td></tr>"
                % (html.escape(a), html.escape(b)) for a, b in rows))
    return head + "".join(parts) + more


def dock(panel):
    """Put the panel INSIDE PyMOL's main window, under the 3D view, so it never
    covers the structure and you can click residues while it is open.  It can
    be dragged to another edge, or pulled out as a separate window, by its
    title bar.  Returns the dock, or None when PyMOL has no Qt main window
    (the panel is then an ordinary window)."""
    try:
        import pymol.gui
        window = pymol.gui.get_qtwindow()
        if window is None or not hasattr(window, "addDockWidget"):
            return None
        holder = QtWidgets.QDockWidget("Prospector for PyMOL", window)
        holder.setObjectName("prospector_panel")
        holder.setAllowedAreas(QtCore.Qt.AllDockWidgetAreas)
        panel.setWindowFlags(QtCore.Qt.Widget)
        scroll = QtWidgets.QScrollArea()          # small screens: scroll, don't squash
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        scroll.setWidget(panel)
        holder.setWidget(scroll)
        window.addDockWidget(QtCore.Qt.BottomDockWidgetArea, holder)
        try:
            # roughly the lower half of the window: big enough to read, and the
            # border above it can be dragged
            window.resizeDocks([holder], [max(430, int(window.height() * 0.48))],
                               QtCore.Qt.Vertical)
        except Exception:
            pass
        return holder
    except Exception:
        return None


class Panel(QtWidgets.QDialog):
    """Ask box on top, the structure you are working on under it, then five
    numbered tabs in the order you would normally work:
        1 Find  ->  2 Understand  ->  3 Look  ->  4 Experiment  ->  5 Keep and learn
    Back / Clear / More and the command that was run stay visible at the bottom."""

    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)
        self.setWindowTitle("Prospector for PyMOL")
        self.resize(1180, 640)
        self.setMinimumSize(560, 300)
        font = self.font()               # a little larger than PyMOL's default
        if font.pointSizeF() > 0:
            font.setPointSizeF(font.pointSizeF() + 1.5)
            self.setFont(font)
        self._row_commands = []     # one PyMOL command per table row

        def field(hint, width=150):
            box = QtWidgets.QLineEdit()
            box.setPlaceholderText(hint)
            box.setMinimumWidth(width)
            return box

        def make_button(text, action, tip=""):
            widget = QtWidgets.QPushButton(text)
            widget.clicked.connect(lambda checked=False: action())
            if tip:
                widget.setToolTip(tip)
            return widget

        def hint(text):
            label = QtWidgets.QLabel(text)
            label.setStyleSheet("color: gray;")
            label.setWordWrap(True)
            # never let a long explanation force the panel wider: it wraps
            label.setSizePolicy(QtWidgets.QSizePolicy.Ignored,
                                QtWidgets.QSizePolicy.Preferred)
            label.setMinimumWidth(120)
            return label

        def page(*rows):
            """A tab page: each row is (title, [widgets], explanation).  The
            explanation sits under its controls, so it has the full width."""
            inner = QtWidgets.QWidget()
            grid = QtWidgets.QGridLayout(inner)
            grid.setColumnStretch(1, 1)
            grid.setVerticalSpacing(2)
            for index, (title, widgets, text) in enumerate(rows):
                name = QtWidgets.QLabel("<b>%s</b>" % title.replace("\n", "<br>"))
                name.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
                grid.addWidget(name, 3 * index, 0, 2, 1)
                line = QtWidgets.QHBoxLayout()
                for item in widgets:
                    line.addWidget(item)
                line.addStretch(1)
                grid.addLayout(line, 3 * index, 1)
                grid.addWidget(hint(text), 3 * index + 1, 1)
                grid.setRowMinimumHeight(3 * index + 2, 8)
            grid.setRowStretch(3 * len(rows), 1)
            area = QtWidgets.QScrollArea()
            area.setWidgetResizable(True)
            area.setFrameShape(QtWidgets.QFrame.NoFrame)
            area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
            area.setWidget(inner)
            return area

        # ---- always on top: ask, and what you are working on -----------------------
        self.ask_box = QtWidgets.QLineEdit()
        self.ask_box.setPlaceholderText("Say what you want - e.g.  load haemoglobin   "
                                        "where does the drug bind?   show it as sticks")
        self.ask_button = make_button("Ask", self.run_ask,
                                      "Free and instant: built-in rules first")
        self.askai_button = make_button("Ask the AI", lambda: self.run_ask(ai=True),
                                        "Straight to the AI (needs a key: prospector ai)")
        ask_row = QtWidgets.QHBoxLayout()
        ask_row.addWidget(QtWidgets.QLabel("<b>Ask</b>"))
        ask_row.addWidget(self.ask_box, 1)
        ask_row.addWidget(self.ask_button)
        ask_row.addWidget(self.askai_button)
        self.holder = None           # the dock this panel sits in, if any
        self._big = False
        self.size_button = make_button("Bigger panel", self.toggle_size,
                                       "Switch between a small and a large panel")
        ask_row.addWidget(self.size_button)

        self.objects = QtWidgets.QComboBox()
        self.objects.setMinimumWidth(160)
        self.objects.setToolTip("The buttons below act on this structure")
        self.working_hint = hint("")
        working = QtWidgets.QHBoxLayout()
        working.addWidget(QtWidgets.QLabel("<b>Working on</b>"))
        working.addWidget(self.objects)
        working.addWidget(self.working_hint, 1)

        # ---- what you clicked in the 3D view, with one-press actions ----------------
        self.click_label = QtWidgets.QLabel("nothing yet - click a residue in the 3D view")
        self.click_explain = make_button("Explain it", lambda: self._click_action("explain"),
                                         "Everything about the clicked residue")
        self.click_zoom = make_button("Zoom to it", lambda: self._click_action("zoom"))
        self.click_change = make_button("Change it...", lambda: self._click_action("change"),
                                        "What if this residue were another amino acid?")
        self.click_measure = make_button("Measure a distance...", self.start_picking,
                                         "Press this, then click two points in the 3D "
                                         "view: the distance between them is measured")
        self.click_forget = make_button("Deselect", lambda: self._click_action("forget"),
                                        "Forget the clicks and remove the pink squares")
        clicked_row = QtWidgets.QHBoxLayout()
        clicked_row.addWidget(QtWidgets.QLabel("<b>Clicked</b>"))
        clicked_row.addWidget(self.click_label)
        for widget in (self.click_explain, self.click_zoom, self.click_change):
            widget.setEnabled(False)
            clicked_row.addWidget(widget)
        clicked_row.addWidget(self.click_measure)
        clicked_row.addWidget(self.click_forget)
        clicked_row.addStretch(1)

        # ---- 1 Find ------------------------------------------------------------------
        self.query = QtWidgets.QLineEdit()
        self.query.setPlaceholderText("a protein or a molecule - e.g. hemoglobin, p53, "
                                      "4HHB, cholesterol, ATP")
        self.find_button = make_button("List matches", self.run_find,
                                       "Show the candidates and why they rank as they do")
        self.grab_button = make_button("Load the best one", self.run_grab,
                                       "Find the best structure and load it")
        self.grab_button.setDefault(True)
        name_row = QtWidgets.QHBoxLayout()
        name_row.addWidget(QtWidgets.QLabel("<b>Name</b>"))
        name_row.addWidget(self.query, 1)
        name_row.addWidget(self.grab_button)
        name_row.addWidget(self.find_button)

        self.organism = QtWidgets.QComboBox()
        self.organism.setEditable(True)
        self.organism.addItems(["", "human", "mouse", "yeast", "e. coli", "sars-cov-2"])
        self.method = QtWidgets.QComboBox()
        self.method.addItems(["", "xray", "em", "nmr"])
        self.resolution = field("e.g. 2.0", 60)
        self.ligand = field("e.g. ATP", 60)
        self.year = field("e.g. 2015", 60)
        for box in (self.resolution, self.ligand, self.year):
            box.setMaximumWidth(80)
        self.filter_box = QtWidgets.QWidget()
        filters = QtWidgets.QHBoxLayout(self.filter_box)
        filters.setContentsMargins(0, 0, 0, 0)
        for label, widget in (("Organism", self.organism), ("Method", self.method),
                              ("Resolution better than", self.resolution),
                              ("Contains ligand", self.ligand),
                              ("Released after", self.year)):
            filters.addWidget(QtWidgets.QLabel(label))
            filters.addWidget(widget)
        filters.addStretch(1)
        self.filter_box.setVisible(False)
        self.filter_toggle = QtWidgets.QCheckBox("Narrow the search (organism, method, "
                                                 "resolution ...)")
        self.filter_toggle.toggled.connect(self.filter_box.setVisible)

        self.table = QtWidgets.QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setMinimumHeight(90)
        self.load_button = make_button("Load the selected row", self.run_selected,
                                       "Or double-click a row")
        self.af_button = make_button("AlphaFold model", lambda: self._run(
            "alphafold %s" % self.arguments()), "The predicted model of this protein")
        self.molecule_button = make_button("As a small molecule", lambda: self._run(
            "molecule %s" % self.query.text().strip()),
            "Load the molecule of this name (cholesterol, ATP ...), not a protein")
        under_table = QtWidgets.QHBoxLayout()
        under_table.addWidget(self.load_button)
        under_table.addWidget(self.af_button)
        under_table.addWidget(self.molecule_button)
        under_table.addStretch(1)
        under_table.addWidget(make_button("Surprise me", lambda: self.run_line("surprise"),
                                          "Load a random classic structure"))
        find_page = QtWidgets.QWidget()
        find_layout = QtWidgets.QVBoxLayout(find_page)
        find_layout.addLayout(name_row)
        find_layout.addWidget(self.filter_toggle)
        find_layout.addWidget(self.filter_box)
        find_layout.addWidget(self.table, 1)
        find_layout.addLayout(under_table)

        # ---- or build one ---------------------------------------------------------
        self.build_kind = QtWidgets.QComboBox()
        for label, kind, hint_text, example in (
                ("DNA double helix", "dna", "the sequence of one strand, 5' to 3' "
                 "(the partner strand is added)", "ATGCGTACGTTA"),
                ("RNA strand", "rna", "the RNA sequence (A, C, G, U)", "GGCUAGCUAC"),
                ("peptide: alpha-helix", "helix", "the sequence in one-letter code",
                 "AKLEELLKKLAE"),
                ("peptide: beta-sheet", "sheet", "two or more strand sequences, "
                 "separated by spaces", "KVTVKVE TVKVTVE"),
                ("peptide: stretched chain", "peptide", "the sequence in one-letter code",
                 "ACDEFGHIK"),
                ("lipid membrane", "membrane", "lipid (popc, dopc, dope), size in nm, "
                 "cholesterol=30 if wanted", "popc 8 cholesterol=30"),
                ("detergent micelle", "micelle", "detergent (sds, ldao) and how many "
                 "molecules", "sds 60"),
                ("whole biological assembly", "assembly", "a protein name or PDB ID - "
                 "leave empty for the structure you are working on", "ferritin")):
            self.build_kind.addItem(label, (kind, hint_text, example))
        self.build_box = field("", 260)
        self.build_hint = hint("")
        build_action = self.run_build
        self.build_box.returnPressed.connect(build_action)
        self.build_kind.currentIndexChanged.connect(lambda *ignored: self._build_hint())
        build_row = QtWidgets.QHBoxLayout()
        build_row.addWidget(QtWidgets.QLabel("<b>Or build a model</b>"))
        build_row.addWidget(self.build_kind)
        build_row.addWidget(self.build_box, 1)
        build_row.addWidget(make_button("Build", build_action,
                                        "Models from rules are marked IDEALISED MODEL; "
                                        "assemblies are real PDB data"))
        find_layout.insertLayout(1, build_row)          # right under the Name row
        find_layout.insertWidget(2, self.build_hint)
        self._build_hint()

        # ---- 2 Understand --------------------------------------------------------------
        self.focus_box = QtWidgets.QComboBox()
        self.focus_box.setEditable(True)
        self.focus_box.setMinimumWidth(230)
        self.focus_box.setToolTip("The key places of this structure. You can also type "
                                  "one: a residue (A/57), a chain (chain B), a molecule")
        self._place_titles = None
        self._asked_places = set()
        self._fill_places()
        self.explain_box = field("click a residue, or type A/57", 190)
        focus_action = lambda: self.run_with("focus %s {name}",
                                             self.focus_box.currentText(),
                                             empty="focus {name}")
        explain_action = lambda: self.run_with("explain {name} %s", self.explain_box.text())
        self.focus_box.activated.connect(lambda *ignored: focus_action())   # pick or Enter
        self.explain_box.returnPressed.connect(explain_action)
        understand_page = page(
            ("What is it?", [make_button("About", lambda: self.run_tool("about"))],
             "The info card: protein, function, method, chains, bound molecules - and "
             "the list of key places."),
            ("Guided tour", [make_button("Start tour", lambda: self.run_tool("tour")),
                             make_button("< Prev", lambda: self.run_line("prev")),
                             make_button("Next >", lambda: self.run_line("next"))],
             "A walk through the structure: fold, domains, functional sites, ligand, "
             "disease variants."),
            ("Zoom to a place", [self.focus_box, make_button("Zoom there", focus_action),
                                 make_button("Whole structure",
                                             lambda: self.run_tool("focus all"))],
             "Pick one of this structure's key places from the list: bound molecules, "
             "sites, domains, chains. The list is filled in a moment after loading."),
            ("One residue", [self.explain_box, make_button("Explain it", explain_action)],
             "Everything about one residue. Tip: click a residue in the 3D window and "
             "its address is filled in for you."),
            ("Known sites", [make_button("Colour the known sites",
                                         lambda: self.run_tool("features"))],
             "Colours what the databases know about this protein: active and binding "
             "sites, domains, disulfide bonds, disease variants."),
            ("Bound molecule", [make_button("Show what holds it",
                                            lambda: self.run_tool("pocket"))],
             "Zooms to the drug or cofactor bound to the protein and shows the residues "
             "and hydrogen bonds that hold it."),
            ("Sequence", [make_button("Show sequence", lambda: self.run_tool("sequence"))],
             "The amino-acid sequence with helices, strands and sites marked."))

        # ---- 3 Look ----------------------------------------------------------------------
        self.styles = QtWidgets.QComboBox()
        for label, value in (("ribbon along the backbone (cartoon)", "cartoon"),
                             ("sticks: every bond", "sticks"),
                             ("balls and sticks", "ballstick"),
                             ("spheres: every atom at full size", "spheres"),
                             ("surface: the outer skin", "surface"),
                             ("thin lines", "lines"),
                             ("thin backbone trace (ribbon)", "ribbon"),
                             ("tube", "tube"), ("mesh surface", "mesh"),
                             ("dotted surface", "dots")):
            self.styles.addItem(label, value)
        self.presets = QtWidgets.QComboBox()
        for label, value in (("each chain its own colour", "chains"),
                             ("rainbow: blue start to red end", "rainbow"),
                             ("helices, strands and loops", "secondary"),
                             ("kind of amino acid", "types"),
                             ("every atom by element", "elements"),
                             ("water-avoiding patches on the surface", "hydrophobicity"),
                             ("flexible parts (thick = moves a lot)", "putty"),
                             ("how buried the bound molecule is", "ligand"),
                             ("clean figure for a report", "publication")):
            self.presets.addItem(label, value)
        self.scout_box = QtWidgets.QComboBox()
        self.scout_box.setEditable(True)
        self.scout_box.addItems(["histidines", "cysteines", "hydrophobic", "charged",
                                 "aromatic", "helices", "strands", "loops", "disulfides",
                                 "salt bridges", "exposed", "buried", "interface",
                                 "flexible", "ligands", "metals"])
        self.choose_box = field("histidines within 5 of the zinc", 250)
        choose_action = lambda: self.run_with("choose %s in {name}", self.choose_box.text())
        self.choose_box.returnPressed.connect(choose_action)
        style_action = lambda: self.run_tool("viewas", after=self.styles.currentData())
        preset_action = lambda: self.run_tool("look %s" % self.presets.currentData())
        scout_action = lambda: self.run_with("scout {name} %s", self.scout_box.currentText())
        # picking from a list does it at once (the button is there to do it again)
        self.styles.activated.connect(lambda *ignored: style_action())
        self.presets.activated.connect(lambda *ignored: preset_action())
        self.scout_box.activated.connect(lambda *ignored: scout_action())
        look_page = page(
            ("Draw it as", [self.styles, make_button("Apply again", style_action)],
             "Pick from the list - the view changes at once. How the atoms are drawn: "
             "ribbon, sticks, balls, surface ..."),
            ("Colour it by", [self.presets, make_button("Apply again", preset_action)],
             "Pick from the list - the colours change at once and the drawing chosen "
             "above is kept: each chain, rainbow, helix/strand, kind of amino acid ..."),
            ("Find all the...", [self.scout_box, make_button("Highlight them",
                                                             scout_action)],
             "Pick from the list - they are highlighted at once, with what they are. "
             "Or type your own (any amino acid works) and press Enter."),
            ("Select some atoms", [self.choose_box, make_button("Select", choose_action)],
             "Describe them in plain words. They are marked with pink squares, and you "
             "are shown the line PyMOL itself would need - a way to learn its syntax."),
            ("Whole view", [make_button("Show everything",
                                        lambda: self.run_tool("focus all")), make_button("White background",
                                                       lambda: self.run_line("bg_color white")),
                            make_button("Black background",
                                        lambda: self.run_line("bg_color black"))],
             "Back to the whole structure, or change the background."))

        # ---- 4 Experiment ------------------------------------------------------------------
        from . import tools as tool_commands
        self.whatif_box = field("residue, e.g. A/248", 110)
        self.whatif_box.setMaximumWidth(130)
        self.whatif_new = QtWidgets.QComboBox()
        for code in sorted(tool_commands.AA, key=lambda c: tool_commands.AA[c][0]):
            name, group, _ = tool_commands.AA[code]
            letter = [k for k, v in tool_commands.ONE_TO_THREE.items() if v == code][0]
            self.whatif_new.addItem("%s  (%s)" % (name, group), letter)
        self.whatif_new.setCurrentIndex(self.whatif_new.findData("A"))
        self.whatif_new.setToolTip("The amino acid it becomes")
        self.whatif_new.setMaximumWidth(210)
        self.measure_a = field("from, e.g. A/176", 110)
        self.measure_b = field("to, e.g. ZN", 110)
        for box in (self.measure_a, self.measure_b):
            box.setMaximumWidth(130)
        self.compare_box = field("another protein or PDB ID, e.g. 1TUP", 240)
        whatif_action = self.run_whatif
        measure_action = self.run_measure
        self.pick_button = make_button("Pick two points by clicking", self.start_picking,
                                       "The next two clicks in the 3D view are measured")
        compare_action = lambda: self.run_with("compare {name}, %s", self.compare_box.text())
        self.whatif_box.returnPressed.connect(whatif_action)
        self.whatif_new.activated.connect(lambda *ignored: self.run_whatif(from_list=True))
        self.measure_b.returnPressed.connect(measure_action)
        self.compare_box.returnPressed.connect(compare_action)
        experiment_page = page(
            ("What if a residue\nwere different?",
             [self.whatif_box, QtWidgets.QLabel("becomes"), self.whatif_new,
              make_button("Show before / after", whatif_action)],
             "1. Click a residue in the 3D view (or type its address).  2. Pick what it "
             "becomes from the list - it runs at once (the button runs it again). You get the original on the LEFT and the "
             "changed protein on the RIGHT, and a list of what changes: charge, size, "
             "contacts, hydrogen bonds. Nothing is changed in the original."),
            ("Back to normal",
             [make_button("One view again", lambda: self.run_line("whatif off"))],
             "After a what-if: leaves the left/right view and shows the original alone."),
            ("How far apart\nare two things?",
             [self.pick_button, QtWidgets.QLabel("  or type:"), self.measure_a,
              QtWidgets.QLabel("to"), self.measure_b,
              make_button("Measure typed", measure_action)],
             "1. Press 'Pick two points by clicking'.  2. Click the first point in the 3D "
             "view, then the second. The distance is measured at once: a dashed yellow "
             "line with the value in Angstrom (under about 4 A they touch). Other clicks "
             "never start a measurement. Or type two addresses (A/176, ZN) and press "
             "'Measure typed'."),
            ("Overlay another\nstructure",
             [self.compare_box, make_button("Compare", compare_action)],
             "Type a second protein. It is loaded and laid on top of this one; red = "
             "parts that differ most, blue = parts that match. Good for normal versus "
             "mutant, or with and without a drug."),
            ("Related proteins",
             [make_button("Find similar proteins", lambda: self.run_tool("similar"))],
             "Lists other proteins with a similar sequence. Nothing changes on screen."))

        # ---- 5 Keep & learn ------------------------------------------------------------------
        from . import teach
        self.session_box = field("a name, e.g. my_work", 150)
        self.quiz_box = field("your guess", 120)
        self.short_box = QtWidgets.QCheckBox("Short answers")
        self.short_box.setChecked(base_commands.short_mode())
        self.short_box.toggled.connect(
            lambda on: self.run_line("prospector short %s" % ("on" if on else "off")))
        self.click_box = QtWidgets.QCheckBox("Clicks fill the boxes")
        self.click_box.setChecked(True)
        self.click_box.setToolTip("Click a residue in the 3D window and its address "
                                  "appears in the Explain, Zoom, Change and Distance boxes")
        self.teach_box = QtWidgets.QCheckBox("Teach mode")
        self.teach_box.setChecked(teach.enabled())
        self.teach_box.toggled.connect(
            lambda on: self.run_line("prospector teach %s" % ("on" if on else "off")))
        keep_page = page(
            ("Picture", [make_button("Save picture...", self.run_snap),
                         make_button("One-page report", lambda: self.run_tool("report"))],
             "A high-resolution PNG of the view, or an HTML page with card and pictures."),
            ("Your work", [self.session_box,
                           make_button("Save", lambda: self.run_line(
                               "session save %s" % (self.session_box.text().strip()
                                                    or "my_work"))),
                           make_button("Open", lambda: self.run_line(
                               "session open %s" % self.session_box.text().strip()
                               if self.session_box.text().strip() else "session"))],
             "Saves everything on screen. Open with an empty box lists what you saved."),
            ("Lists", [make_button("History", lambda: self.run_line("history")),
                       make_button("Add favourite", lambda: self.run_tool("fav add")),
                       make_button("Favourites", lambda: self.run_line("fav")),
                       make_button("Load a list...", self.run_batch)],
             "Past searches, your favourites, or load many structures from a text file."),
            ("Quiz", [make_button("New mystery structure", lambda: self.run_line("quiz")),
                      self.quiz_box,
                      make_button("Guess", lambda: self.run_line(
                          "quiz guess %s" % self.quiz_box.text().strip())
                          if self.quiz_box.text().strip() else self.message.setText(
                              "Type your guess in the box first, e.g. myoglobin.")),
                      make_button("Hint", lambda: self.run_line("quiz hint")),
                      make_button("Give up", lambda: self.run_line("quiz reveal"))],
             "A structure is loaded without its name. Look at it, type what you think "
             "it is, and press Guess."),
            ("Settings", [self.short_box, self.teach_box, self.click_box],
             "Short answers: PyMOL's text area gets only the answer. Teach mode: also "
             "prints the PyMOL commands behind every action, so you learn them."),
            ("Help and AI", [make_button("List all commands", lambda: self.run_line("info")),
                             make_button("First steps", lambda: self.run_line("start")),
                             make_button("AI: status", lambda: self.run_line("prospector ai")),
                             make_button("AI: enter my key",
                                         lambda: self.run_line("prospector ai key"))],
             "The AI is optional: it reads sentences the Ask box does not understand. "
             "It needs your own key."))

        self.tabs = QtWidgets.QTabWidget()
        for title, widget, tip in (
                ("1  Find", find_page, "Find and load a protein or a molecule"),
                ("2  Understand", understand_page, "What it is and where its parts are"),
                ("3  Look", look_page, "How it is drawn and coloured; spot and select"),
                ("4  Experiment", experiment_page, "Change, measure, compare"),
                ("5  Keep and learn", keep_page, "Save, practise, settings")):
            self.tabs.addTab(widget, title)
            self.tabs.setTabToolTip(self.tabs.count() - 1, tip)

        # ---- always at the bottom ----------------------------------------------------------
        self.preview = QtWidgets.QLabel()
        self.preview.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.preview.setToolTip("Every button runs a command you could type yourself")
        self.message = QtWidgets.QLabel("Start in tab 1: type a name and press Enter. Or "
                                        "just ask in the box at the top. The panel can "
                                        "stay open: click residues in the 3D view above.")
        self.message.setWordWrap(True)
        self.message.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        bottom = QtWidgets.QHBoxLayout()
        bottom.addWidget(self.preview, 1)
        for text, line, tip in (("Undo last step", "back",
                                 "Go back one step: the picture, and anything the step "
                                 "added"),
                                ("Clear everything", "clear",
                                 "Empty the 3D window ('Undo last step' right after "
                                 "brings it back)")):
            bottom.addWidget(make_button(text, lambda l=line: self.run_line(l), tip))

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(ask_row)
        layout.addLayout(working)
        layout.addLayout(clicked_row)
        self.result = QtWidgets.QTextEdit()
        self.result.setReadOnly(True)
        self.result.setHtml("<b>WHAT JUST HAPPENED</b><p>Press a button and its result "
                            "appears here in plain words: what is now on screen and what "
                            "it means.</p><p style='color:gray'>The same text is printed "
                            "in PyMOL's own text area.</p>")
        self._last_output = []
        self.more_button = make_button("Show all of it", self.show_full_output,
                                       "The complete text of the last result")
        self.more_button.setEnabled(False)
        result_side = QtWidgets.QWidget()
        result_layout = QtWidgets.QVBoxLayout(result_side)
        result_layout.setContentsMargins(0, 0, 0, 0)
        result_layout.addWidget(self.result, 1)
        result_layout.addWidget(self.more_button)
        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.splitter.addWidget(self.tabs)
        self.splitter.addWidget(result_side)
        self.splitter.setChildrenCollapsible(False)
        self.tabs.setMinimumWidth(420)
        result_side.setMinimumWidth(320)
        self.splitter.setStretchFactor(0, 3)
        self.splitter.setStretchFactor(1, 2)
        self.splitter.setSizes([760, 440])
        layout.addWidget(self.splitter, 1)
        layout.addWidget(self.message)
        layout.addLayout(bottom)

        self.ask_box.returnPressed.connect(self.run_ask)
        self.query.returnPressed.connect(self.run_grab)
        self.table.cellDoubleClicked.connect(lambda row, column: self.run_selected())
        for signal in (self.query.textChanged, self.organism.editTextChanged,
                       self.method.currentIndexChanged, self.resolution.textChanged,
                       self.ligand.textChanged, self.year.textChanged):
            signal.connect(lambda *ignored: self._update_preview())
        self.refresh_objects()
        self._update_preview()
        # keep "Working on" in step with PyMOL, whatever loaded or removed things
        self._seen = None
        self._clicks = []            # the clicked residues, newest last
        self._picking = False        # True while the next two clicks are being measured
        self._picked = []
        self._selection = self._selected_residues()
        self._timer = QtCore.QTimer(self)
        self._timer.timeout.connect(self._keep_in_step)
        self._timer.start(500)

    def _keep_in_step(self):
        try:
            names = tuple(n for n in (cmd.get_names("objects") or [])
                          if cmd.get_type(n) == "object:molecule")
        except Exception:
            return
        if names != self._seen:
            self._seen = names
            self.refresh_objects()
        self._fill_places()
        self._take_click()

    def _fill_places(self):
        """The 'Zoom to a place' list = the key places of the structure worked on."""
        from . import places
        obj = self.objects.currentText() if hasattr(self, "objects") else ""
        listed = places._places.get(obj) or []
        if obj and not listed and obj not in self._asked_places:
            # work out this structure's key places quietly, in the background
            self._asked_places.add(obj)
            identity = base_commands._identify(obj)
            if identity and identity[0] != "mol":
                def built(result, error, name=obj):
                    try:
                        if not error and result and name in base_commands._molecules():
                            places.build(name, result[0], result[1])
                    except Exception:
                        pass
                from . import worker
                worker.submit(lambda: places.fetch_annotation(identity), built)
        titles = tuple(place.title for place in listed) or (
            "ligand", "binding sites", "active site", "chain A", "all")
        if (obj, titles) == self._place_titles:
            return
        self._place_titles = (obj, titles)
        typed = self.focus_box.currentText()
        self.focus_box.blockSignals(True)
        self.focus_box.clear()
        self.focus_box.addItems(list(titles))
        if typed and typed in titles:
            self.focus_box.setCurrentIndex(titles.index(typed))
        self.focus_box.blockSignals(False)

    def _selected_residues(self):
        """Every residue in PyMOL's click selection "sele", as a set."""
        found = set()
        try:
            if "sele" in cmd.get_names("selections"):
                cmd.iterate("sele", "found.add((model, chain, resi, resn))",
                            space={"found": found})
        except Exception:
            pass
        return found

    def _take_click(self):
        """Notice a residue clicked in the 3D window.  PyMOL ADDS each click to
        its selection, so the click is the residue that is new since last time."""
        if not self.click_box.isChecked():
            return
        now = self._selected_residues()
        new = now - self._selection
        self._selection = now
        if len(new) != 1:
            return                     # nothing new, a deselect, or a whole region
        obj, chain, resi, resn = list(new)[0]
        polymer = False
        try:
            polymer = cmd.count_atoms("sele and polymer and model %s and resi %s%s" % (
                obj, ("\\" + resi) if resi.startswith("-") else resi,
                " and chain %s" % chain if chain else "")) > 0
        except Exception:
            pass
        spec = "%s/%s" % (chain, resi) if chain else resi
        end = spec if polymer else resn            # a bound molecule goes by its code
        base_commands._last_click[0] = (obj, chain, resi, resn)
        index = self.objects.findText(obj)
        if index >= 0:
            self.objects.setCurrentIndex(index)
        self._clicks = (self._clicks + [{"end": end, "spec": spec, "polymer": polymer,
                                         "name": "%s %s" % (resn, resi), "obj": obj}])[-2:]
        self._show_clicks()

    def _show_clicks(self):
        """What a click does.  Normally it names the residue and fills the
        Explain / Zoom / Change boxes.  While measuring (after "Measure a
        distance" was pressed) the next two clicks are the two points."""
        last = self._clicks[-1]
        if self._picking:
            self._picked.append(last)
            if len(self._picked) == 1:
                self.measure_a.setText(last["end"])
                self.measure_b.setText("")
                self._measure_hint("Point 1 = %s. Now click point 2." % last["name"])
                return
            first, second = self._picked[0], self._picked[1]
            self.measure_b.setText(second["end"])
            self._stop_picking()
            if first["obj"] != second["obj"]:
                self.message.setText("Those two are in different structures (%s and %s); "
                                     "a distance needs both in the same one. Press the "
                                     "measure button to try again."
                                     % (first["obj"], second["obj"]))
                return
            self.click_label.setText("<b>measured %s - %s</b>" % (first["name"],
                                                                   second["name"]))
            self.run_line("measure %s %s %s" % (first["obj"], first["end"], second["end"]))
            return
        if last["polymer"]:
            self.explain_box.setText(last["spec"])
            self.focus_box.setEditText(last["spec"])
            self.whatif_box.setText(last["spec"])
        for button in (self.click_explain, self.click_zoom, self.click_change):
            button.setEnabled(last["polymer"])
        self.click_label.setText("<b>%s</b>" % last["name"])
        self.message.setText("You clicked %s in %s. Use the buttons in the Clicked row. "
                             "(To measure, press 'Measure a distance' first.)"
                             % (last["name"], last["obj"]))

    # ---- measuring: press the button, then click two points -----------------
    def start_picking(self):
        """The next two clicks in the 3D view are the two ends of a distance."""
        if self._picking:                       # pressed again = cancel
            self._stop_picking()
            self.click_label.setText("measuring cancelled - click a residue in the 3D view")
            self.message.setText("Measuring cancelled.")
            return
        if not self.click_box.isChecked():
            self.click_box.setChecked(True)
        self._picking = True
        self._picked = []
        self.measure_a.setText("")
        self.measure_b.setText("")
        try:                                    # start from a clean selection
            cmd.delete("sele")
        except Exception:
            pass
        self._selection = set()
        for button in (self.click_measure, self.pick_button):
            button.setText("Cancel measuring")
        self._measure_hint("Measuring: click point 1 in the 3D view.")

    def _stop_picking(self):
        self._picking = False
        self._picked = []
        self.click_measure.setText("Measure a distance...")
        self.pick_button.setText("Pick two points by clicking")
        self.click_label.setStyleSheet("")

    def _measure_hint(self, text):
        self.click_label.setText("<b>%s</b>" % text)
        self.click_label.setStyleSheet("color: #b35900;")
        self.message.setText(text + "  (Press 'Cancel measuring' to stop.)")

    def _click_action(self, what):
        if not self._clicks and what != "forget":
            return
        if not self._clicks:
            self._clicks = [{"obj": "", "spec": "", "name": "", "end": "", "polymer": False}]
        last = self._clicks[-1]
        if what == "explain":
            self.run_line("explain %s %s" % (last["obj"], last["spec"]))
        elif what == "zoom":
            self.run_line("focus %s %s" % (last["spec"], last["obj"]))
        elif what == "change":
            self.tabs.setCurrentIndex(3)
            self.whatif_box.setText(last["spec"])
            self.whatif_box.setFocus()
            self.whatif_box.setCursorPosition(len(last["spec"]))
            self.message.setText("Now pick what %s becomes in the list next to it, then "
                                 "press 'Show before / after'." % last["name"])
        elif what == "forget":
            if self._picking:
                self._stop_picking()
            self._clicks = []
            self._selection = set()
            try:
                cmd.deselect()
                cmd.delete("sele")
            except Exception:
                pass
            self.click_label.setText("nothing yet - click a residue in the 3D view")
            for button in (self.click_explain, self.click_zoom, self.click_change):
                button.setEnabled(False)

    # ---- building the command line -------------------------------------
    def arguments(self):
        parts = [self.query.text().strip()]
        organism = self.organism.currentText().strip()
        if organism:
            parts.append('organism="%s"' % organism if " " in organism
                         else "organism=%s" % organism)
        if self.method.currentText():
            parts.append("method=%s" % self.method.currentText())
        if self.resolution.text().strip():
            parts.append("res<%s" % self.resolution.text().strip())
        if self.ligand.text().strip():
            parts.append("ligand=%s" % self.ligand.text().strip())
        if self.year.text().strip():
            parts.append("year>%s" % self.year.text().strip())
        return " ".join(p for p in parts if p)

    def _update_preview(self):
        self.preview.setText("Command:  grab %s" % self.arguments())

    def _run(self, command):
        if not self.query.text().strip():
            self.message.setText("Type a protein name first.")
            return
        self.message.setText("Searching... (PyMOL stays usable meanwhile)")
        self.preview.setText("Command:  %s" % command)
        cmd.do(command)

    def run_ask(self, ai=False):
        text = self.ask_box.text().strip()
        if not text:
            self.message.setText("Type what you want in the box at the top, e.g.  show "
                                 "the active site")
            return
        command = "%s %s" % ("askai" if ai else "ask", text)
        self.preview.setText("Command:  %s" % command)
        self.message.setText("See the PyMOL console: it shows which command your "
                             "sentence became, then runs it.")
        cmd.do(command)

    def run_find(self):
        self._run("find %s" % self.arguments())

    def run_grab(self):
        self._run("grab %s" % self.arguments())

    def refresh_objects(self):
        """Fill the drop-down with the molecules currently loaded in PyMOL."""
        current = self.objects.currentText()
        names = [n for n in (cmd.get_names("objects") or [])
                 if cmd.get_type(n) == "object:molecule"]
        self.objects.clear()
        self.objects.addItems(names)
        if current in names:
            self.objects.setCurrentIndex(names.index(current))
        elif names:
            self.objects.setCurrentIndex(len(names) - 1)
        self.working_hint.setText(
            "the buttons in tabs 2-5 act on this structure" if names else
            "nothing loaded yet - type a name in tab 1; tabs 2-5 will load it for you")

    def run_line(self, line):
        self.preview.setText("Command:  %s" % line)
        self.message.setText("Running:  %s   - the result appears in the box on the "
                             "right." % line)
        cmd.do(line)

    # ---- the result box, and the panel's size ----------------------------
    def on_output(self, command, lines):
        """Called when any Prospector command has finished printing."""
        lines = [line for line in lines if line.strip()]
        if not lines or command.split()[:1] == ["more"]:
            return
        self._last_output = ["You ran:  %s" % command, ""] + lines
        kept, hidden = base_commands._shorten(
            lines, blocks_kept=3, lines_each=5,
            always=("IN SHORT", "WHAT IT MEANS", "WHAT IT IS", "WHEN DONE", "CHOSEN",
                    "PYMOL"))
        if hidden <= 3:
            kept, hidden = lines, 0
        self._last_command = command
        self.result.setHtml(_readable(command, kept, hidden))
        if command.split()[:1] and command.split()[0] in ("build", "grab", "molecule",
                                                          "alphafold", "surprise"):
            self.refresh_objects()                # follow what was just made or loaded
            index = self.objects.findText(base_commands._last[0] or "")
            if index >= 0:
                self.objects.setCurrentIndex(index)
        self.more_button.setEnabled(bool(hidden))

    def show_full_output(self):
        if self._last_output:
            self.result.setHtml(_readable(getattr(self, "_last_command", ""),
                                          self._last_output[2:], 0))
            self.more_button.setEnabled(False)

    def toggle_size(self):
        self.set_size(not self._big)

    def set_size(self, big):
        """Small: a strip under the 3D view.  Big: most of the window."""
        self._big = bool(big)
        self.size_button.setText("Smaller panel" if self._big else "Bigger panel")
        try:
            if self.holder is not None and not self.holder.isFloating():
                window = self.holder.parent()
                share = 0.68 if self._big else 0.40
                window.resizeDocks([self.holder], [int(window.height() * share)],
                                   QtCore.Qt.Vertical)
            else:
                target = self.holder if self.holder is not None else self
                target.resize(1250 if self._big else 980, 820 if self._big else 520)
        except Exception:
            pass

    def _build_hint(self):
        kind, hint_text, example = self.build_kind.currentData()
        self.build_box.setPlaceholderText("e.g.  %s" % example)
        self.build_hint.setText("Type %s, then press Build. %s" % (
            hint_text, "Nothing to type: Build makes the default." if kind in
            ("membrane", "micelle") else ""))

    def run_build(self):
        kind, _, example = self.build_kind.currentData()
        text = " ".join(self.build_box.text().split())
        if not text and kind in ("dna", "rna", "helix", "sheet", "peptide"):
            self.message.setText("Type a sequence in the box first, e.g. %s" % example)
            return
        if not text and kind == "assembly":
            self.refresh_objects()
            text = self.objects.currentText()
            if not text:
                self.message.setText("Type a protein name or PDB ID, e.g. ferritin")
                return
        self.run_line(" ".join(("build %s %s" % (kind, text)).split()))

    def run_whatif(self, from_list=False):
        """Build the change from the residue box and the "becomes" list."""
        spec = "".join(self.whatif_box.text().split())
        if not spec:
            self.message.setText("First say which residue: click one in the 3D view, or "
                                 "type its address (e.g. A/248) in the box. Then pick "
                                 "what it becomes - it runs at once.")
            self.whatif_box.setFocus()
            return
        if from_list and spec[-1].isalpha() and spec[-2:-1].isdigit():
            spec = spec[:-1]                  # R248Q typed, then another pick: R248 + new
        if not spec[-1].isalpha() or not any(ch.isdigit() for ch in spec):
            spec += self.whatif_new.currentData()         # A/248 + W
        elif spec[-1].isalpha() and spec[-2:-1].isdigit():
            pass                                           # already complete: R248Q
        self.run_with("whatif {name} %s", spec)

    def run_measure(self):
        first, second = self.measure_a.text().strip(), self.measure_b.text().strip()
        if not first or not second:
            self.message.setText("Type both ends first (e.g. A/176 and ZN) - or press "
                                 "'Pick two points by clicking' and click them.")
            return
        self.run_with("measure {name} %s", "%s %s" % (first, second))

    def _target(self):
        """What the tool buttons act on: the loaded structure chosen at the top,
        or - when nothing is loaded - the name typed in tab 1 (the command then
        finds and loads it by itself)."""
        self.refresh_objects()
        name = self.objects.currentText() or self.query.text().strip()
        if not name:
            self.tabs.setCurrentIndex(0)
            self.query.setFocus()
            self.message.setText("Nothing is loaded yet. Type a protein or molecule name "
                                 "here in tab 1 - any button will then load it for you.")
        return name

    def run_tool(self, command, name_only_id=False, after=""):
        name = self._target()
        if name:
            self.run_line(" ".join(("%s %s %s" % (command, name, after)).split()))

    def run_with(self, template, text, empty=None):
        """Run a command that needs the loaded structure and something typed in
        a box, e.g. template "scout {name} %s" with text "histidines"."""
        name = self._target()
        if not name:
            return
        text = " ".join(text.split())
        if not text:
            if empty:
                self.run_line(empty.replace("{name}", name))
            else:
                self.message.setText("Fill in the box next to that button first.")
            return
        self.run_line(" ".join((template % text).replace("{name}", name).split()))

    def run_snap(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save picture as", "",
                                                        "PNG image (*.png)")
        if path:
            self.run_line("snap %s" % path)

    def run_batch(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose a list of proteins",
                                                        "", "Text files (*.txt);;All (*)")
        if path:
            self.run_line("batch %s" % path)

    def run_selected(self):
        row = self.table.currentRow()
        if 0 <= row < len(self._row_commands):
            command = self._row_commands[row]
            self.preview.setText("Command:  %s" % command)
            self.message.setText("Running:  %s" % command)
            cmd.do(command)

    # ---- showing results ------------------------------------------------
    def _fill(self, rows, commands):
        self.table.setRowCount(len(rows))
        for r, row in enumerate(rows):
            for c, value in enumerate(row):
                self.table.setItem(r, c, QtWidgets.QTableWidgetItem(str(value)))
        self._row_commands = commands
        self.table.resizeColumnsToContents()
        if rows:
            self.table.selectRow(0)

    def on_result(self, command_name, result):
        """Called by the command layer whenever find/grab finishes."""
        self.refresh_objects()
        if result.kind == "molecule":
            molecule = result.molecule
            self._fill([[molecule.comp_id or "PubChem", "", "molecule", "-", "", "", "", "",
                         "%s - %s, %s" % (molecule.title(), molecule.kind(),
                                          molecule.formula)]],
                       ["molecule %s" % (molecule.comp_id or molecule.name)])
            self.message.setText("'%s' is a small molecule. Double-click to load it; for "
                                 "proteins named after it add  kind=protein ."
                                 % result.query.text)
            return
        if result.kind == "results":
            rows, commands = [], []
            for hit in result.hits:
                entry = hit.entry
                rows.append([
                    entry.pdb_id + (" *" if hit.classic else ""), "%.1f" % hit.score,
                    entry.method_short(),
                    "%.2f" % entry.resolution if entry.resolution is not None else "-",
                    entry.year or "?", entry.organism_text(hit.accession),
                    entry.chain_count or "?", ",".join(entry.ligands[:4]) or "-", entry.title])
                commands.append("grab %s" % entry.pdb_id)
            self._fill(rows, commands)
            if command_name == "grab":
                self.message.setText("Loaded %s. Next: tab 2 'Understand' tells you what "
                                     "it is - or click any residue in the 3D view."
                                     % result.top.pdb_id)
            else:
                self.message.setText("%d matches, best first. Double-click a row to load "
                                     "it. %s" % (len(result.hits), rank.explain(result.hits)
                                                 if result.mode != "id" else ""))
        elif result.kind == "ambiguous":
            rows, commands = [], []
            choices = result.resolution.choices
            for group in choices:
                command = "find %s" % " ".join([group.suggested_query(choices)]
                                               + result.query.filter_words())
                rows.append([group.lead.gene or group.lead.accession, "", "", "", "",
                             group.organisms_text(), "", "",
                             "%s - %d structures" % (group.label(), group.count)])
                commands.append(command)
            self._fill(rows, commands)
            self.message.setText("'%s' could mean several proteins. Double-click the one "
                                 "you want." % result.resolution.name)
        elif result.kind == "suggest":
            self._fill([["?", "", "", "", "", "", "", "", "Did you mean '%s'? (double-click)"
                         % result.suggestion]], [result.command()])
            self.message.setText(result.message)
        elif result.kind == "alphafold":
            model = result.model
            self._fill([[model.model_id, "", "AlphaFold", "-", "", model.organism, 1, "-",
                         "Predicted model of %s (no experimental structure)" % model.name]],
                       [result.command()])
            self.message.setText("No experimental structure exists; an AlphaFold "
                                 "prediction is available. Double-click to load it.")
        else:
            self._fill([], [])
            self.message.setText(result.message or "Nothing found.")
