"""Command-line entry point: ``pyprimer serve``, quick primer analysis and AutoDimer."""

from __future__ import annotations

import argparse
import sys
import webbrowser
from collections.abc import Sequence

from . import __version__
from .autodimer import AutoDimerSettings, parse_oligos, report_text, screen_dimers, screen_hairpins
from .report import analyse_pair


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pyprimer",
        description="Primer design for PCR, bisulphite PCR, real-time PCR and sequencing.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    serve = commands.add_parser("serve", help="start the web interface")
    serve.add_argument(
        "--host", default="127.0.0.1", help="interface to bind (default: %(default)s)"
    )
    serve.add_argument("--port", type=int, default=8000, help="port (default: %(default)s)")
    serve.add_argument("--open", action="store_true", help="open the interface in a browser")
    serve.add_argument("--debug", action="store_true", help="run Flask's debug server")

    tm = commands.add_parser("analyse", help="print Tm and primer-dimers for one or two primers")
    tm.add_argument("forward")
    tm.add_argument("reverse", nargs="?")

    defaults = AutoDimerSettings()
    ad = commands.add_parser(
        "autodimer", help="screen a set of oligos for primer-dimers or hairpins (AutoDimer)"
    )
    ad.add_argument("file", help="text/FASTA file of oligos ('-' for standard input)")
    ad.add_argument("--hairpin", action="store_true", help="screen for hairpins instead")
    ad.add_argument(
        "--score", type=int, default=defaults.min_score, help="minimum score (default: %(default)s)"
    )
    ad.add_argument(
        "--na", type=float, default=defaults.sodium_m, help="Na+ in M (default: %(default)s)"
    )
    ad.add_argument(
        "--ct",
        type=float,
        default=defaults.strand_um,
        help="total strand concentration in uM (default: %(default)s)",
    )
    ad.add_argument(
        "--dg-temp",
        type=float,
        default=defaults.dg_temp_c,
        help="temperature for dG in C (default: %(default)s)",
    )

    args = parser.parse_args(argv)
    if args.command == "serve":
        return _serve(args.host, args.port, args.open, args.debug)
    if args.command == "autodimer":
        return _autodimer(args)
    return _analyse(args.forward, args.reverse)


def _serve(host: str, port: int, open_browser: bool, debug: bool) -> int:
    from .web import create_app

    app = create_app()
    url = f"http://{'localhost' if host in ('127.0.0.1', '0.0.0.0') else host}:{port}/"
    print(f"PyPrimer {__version__} running at {url} (press Ctrl+C to stop)")
    if open_browser:
        webbrowser.open(url)
    if debug:
        app.run(host=host, port=port, debug=True)
        return 0
    try:
        from waitress import serve
    except ImportError:
        app.run(host=host, port=port, threaded=True)
    else:
        serve(app, host=host, port=port)
    return 0


def _analyse(forward: str, reverse: str | None) -> int:
    analysis = analyse_pair(forward, reverse)
    for error in analysis.errors:
        print(f"Error: {error}", file=sys.stderr)
    for name, primer in (("Forward", analysis.forward), ("Reverse", analysis.reverse)):
        if primer:
            print(
                f"{name}: {primer.seq}\n  Tm {primer.tm:.2f} °C, {primer.length} bases, "
                f"{primer.gc:.0f}% GC, dH {primer.dh:.2f} kcal/mol, dS {primer.ds:.2f} eu, "
                f"dG {primer.dg:.2f} kcal/mol"
            )
    print()
    print(analysis.dimer_text())
    return 1 if analysis.errors else 0


def _autodimer(args: argparse.Namespace) -> int:
    try:
        if args.file == "-":
            text = sys.stdin.read()
        else:
            with open(args.file, encoding="utf-8", errors="replace") as handle:
                text = handle.read()
        settings = AutoDimerSettings(args.na, args.ct, args.dg_temp, args.score)
        oligos = parse_oligos(text)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    if len(oligos) < (1 if args.hairpin else 2):
        print("Error: not enough sequences to screen", file=sys.stderr)
        return 1
    screen = (screen_hairpins if args.hairpin else screen_dimers)(oligos, settings)
    print(report_text(screen), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
