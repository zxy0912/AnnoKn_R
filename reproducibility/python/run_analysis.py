"""Python adapter for run.R; uses its resolved configuration and preflight."""
import ast
import csv
from pathlib import Path, PureWindowsPath
import runpy
import sys


def references(script):
    for node in ast.walk(ast.parse(Path(script).read_text(), filename=script)):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "repro_path":
            args = [ast.literal_eval(arg) for arg in node.args]
            if not 1 <= len(args) <= 2 or not all(isinstance(arg, str) for arg in args):
                raise ValueError("Path keys and suffixes must be literal strings")
            yield args[0], args[1] if len(args) > 1 else ""


def main():
    if sys.argv[1] == "--references":
        for key, relative in references(sys.argv[2]):
            print(key + "\t" + relative)
        return
    if sys.argv[1] != "--paths":
        raise SystemExit("Use Rscript reproducibility/run.R relative/script.py [arguments]")
    with open(sys.argv[2], newline="") as handle:
        paths = {row["key"]: row["path"] for row in csv.DictReader(handle, delimiter="\t")}

    def repro_path(key, relative=None):
        root = paths[key]
        if relative is None:
            return root
        if (not isinstance(relative, str) or Path(relative).is_absolute()
                or PureWindowsPath(relative).drive or ".." in relative.replace("\\", "/").split("/")):
            raise ValueError("Suffix must be relative without '..'")
        return root + "/" + relative

    script = sys.argv[3]
    sys.argv = [script] + sys.argv[4:]
    runpy.run_path(script, init_globals={"repro_path": repro_path}, run_name="__main__")


if __name__ == "__main__":
    main()
