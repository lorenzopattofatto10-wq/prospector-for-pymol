"""Build dist/prospector.zip - the file you give to PyMOL's Plugin Manager.

    python build_zip.py
"""
import os
import zipfile

ROOT = os.path.dirname(os.path.abspath(__file__))


def main():
    os.makedirs(os.path.join(ROOT, "dist"), exist_ok=True)
    target = os.path.join(ROOT, "dist", "prospector.zip")
    count = 0
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
        for folder, subfolders, files in os.walk(os.path.join(ROOT, "prospector")):
            subfolders[:] = [s for s in subfolders if s != "__pycache__"]
            for name in sorted(files):
                if name.endswith((".pyc", ".pyo")):
                    continue
                path = os.path.join(folder, name)
                archive.write(path, os.path.relpath(path, ROOT))
                count += 1
        for extra in ("LICENSE", "LICENSE-docs.md", "README.md", "CHANGELOG.md"):   # credits travel with it
            path = os.path.join(ROOT, extra)
            if os.path.exists(path):
                archive.write(path, os.path.join("prospector", extra))
                count += 1
    print("Wrote %s (%d files, %.0f kB)" % (target, count, os.path.getsize(target) / 1000.0))


if __name__ == "__main__":
    main()
