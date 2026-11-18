# Contributing to Prospector for PyMOL

Thank you for helping. Prospector is meant to make protein structures easy to
explore for students, so reports from people who are *new* to PyMOL are
especially valuable.

## Report a problem or ask for help

Open an issue: https://github.com/lorenzopattofatto10-wq/prospector-for-pymol/issues

Please include:

- your PyMOL version (shown when PyMOL starts) and operating system;
- the Prospector version (the log line `Prospector x.y.z loaded`);
- the exact command you typed and what Prospector printed – copy the text from
  the PyMOL console rather than describing it;
- what you expected instead.

Questions ("how do I…?") are welcome as issues too.

## Suggest an improvement

Open an issue describing what you wanted to do and what got in the way. A real
example ("I typed `ask colour the active site red` and nothing happened") is
more useful than a general wish.

## Contribute code or teaching material

1. Fork the repository and create a branch.
2. Make your change. Keep the two layers apart: `prospector/core/` is pure
   Python (no PyMOL), `prospector/pymol_layer/` holds the PyMOL commands and
   the panel.
3. Run the tests – they work offline:

   ```bash
   python -m unittest discover -s tests -v
   ```

   If you add or change behaviour, add a test for it.
4. Open a pull request explaining what changed and why.

Teaching materials (worksheets like those in `docs/worksheets/`) are very
welcome: please test every command in PyMOL before submitting.

## Ground rules

- Be kind and patient – many users are learning both biology and software.
- Never put API keys or personal data in issues, code or tests.
- Contributions are released under the project's [MIT licence](LICENSE).
