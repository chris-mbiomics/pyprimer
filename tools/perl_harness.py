"""Build a headless Perl harness from the original PerlPrimer source.

The harness contains PerlPrimer's own algorithm subroutines (copied verbatim
from ``perlprimer.pl``) with the Tk GUI stubbed out. It reads JSON requests
on stdin, one per line, and writes one JSON response per line, so the Python
port can be compared against the original implementation.

Usage::

    python tools/perl_harness.py path/to/perlprimer.pl > harness.pl
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

# Subroutines whose bodies are copied from perlprimer.pl.
SUBS = [
    "complement", "bisulphite", "primer_window", "check_degenerate", "bisul_window",
    "calc_amplicon", "calc_seq_primers", "gc", "cc", "tm", "recalculate_dG",
    "primer_dimer", "draw_dimer", "check_range", "find_orf", "find_cpg", "load_data",
]  # fmt: skip

PRELUDE = r"""
use strict;
no warnings;
use JSON::PP;
$| = 1;

our (%oligo_dH, %oligo_dH_full, %oligo_dS, %oligo_dS_full, %genetic_code);
our ($oligo_conc, $mg_conc, $monovalent_cation_conc, $dntp_conc) = (200, 1.5, 50, 0.2);
our ($pd_extensible, $pd_temperature) = (1, 37);
our ($exclude_gc, $exclude_clamp, $min_gc, $max_gc) = (1, 1, 40, 60);
our ($exclude_rr, $repeat, $run) = (1, 3, 4);
our ($exclude_rr_bs, $repeat_bs, $run_bs) = (1, 3, 4);
our ($bisul_min_c, $pre_bs, $exclude_3c, $exclude_cpg) = (30, 1, 1, 0);
our ($ie_overlap, $exclude_ie, $ie_span) = (1, 7, 1);
our ($ie_limit, $ie_limit_5p, $ie_limit_3p);
our ($seq_spacing_min, $seq_spacing_max, $exclude_pd_seq, $seq_pd_min) = (500, 700, 1, 5);
our ($cpg_window, $min_cpg_island, $cpg_gc, $cpg_oe, $cpgplot_method) = (200, 200, 50, 0.6, 0);
our ($bs, $qpcr_flag, $cancel, $reverse, $pfkeys, $pfl, $prl, $pl, $rprimer_r) = (0, 0, 0);
our (@PF, @PR, @primer_pairs, @intron_exon_bounds, @score, @score_sort, @bind_string);
our (%primer_hash, %rating_hash, %packed_widgets, %vars);
our %oligo_dG = (initC => 0.98, initG => 0.98, initA => 1.03, initT => 1.03);

package Stub { sub new { bless {}, shift } sub AUTOLOAD { } sub DESTROY { } }
$packed_widgets{bisul_seq} = Stub->new;
our $seq_widget = Stub->new;

sub sbarprint { }
sub dialogue { push @{ $vars{messages} }, $_[0] }
sub which_nb_page { 'pd' }
sub clean_seq { my $s = shift; $s =~ s/[\s\d-]//g; return $s }
sub get_variables {
    return map { $_ eq 'seq' ? \$seq_widget : (exists $vars{$_} ? \$vars{$_} : \do { my $x }) } @_;
}
"""

DRIVER = r"""
load_data();
recalculate_dG();

sub pairs_out {
    my @sorted = sort { $b->[13] <=> $a->[13] } @primer_pairs;
    @sorted = sort { $b->[10] <=> $a->[10] } @sorted;
    return [ map { [ @$_[0 .. 10, 13] ] } @sorted ];
}

sub design {
    my ($r, $window_sub, $bisulphite) = @_;
    %vars = (%vars, %{ $r->{vars} });
    $vars{messages} = [];
    my $seq = uc $r->{seq};
    check_range();
    $vars{max_ampsize} = $vars{max_range_3p} - $vars{max_range_5p}
        if $vars{max_range_5p} && $vars{max_range_3p};
    if ($vars{min_range} && $vars{max_range}) {
        $vars{min_ampsize} = $vars{max_range} - $vars{min_range};
        $vars{max_ampsize} = $vars{max_range} - $vars{min_range}
            if $vars{max_range} - $vars{min_range} > $vars{max_ampsize};
    }
    my $rev = reverse $seq;
    $reverse = 0;
    @PF = $window_sub->($seq, length $seq);
    $pfkeys = @PF;
    $reverse = 1;
    @PR = $window_sub->($bisulphite ? $rev : complement($rev), length $seq);
    calc_amplicon(length $seq);
    return pairs_out();
}

my $json = JSON::PP->new->canonical;
while (my $line = <STDIN>) {
    my $r = $json->decode($line);
    my $op = $r->{op};
    %vars = ();
    $qpcr_flag = 0;
    @intron_exon_bounds = ();
    my $out;
    if ($op eq 'tm') {
        $out = [ tm($r->{seq}) ];
    } elsif ($op eq 'dimer') {
        my $dg = primer_dimer($r->{f}, $r->{r}, $r->{full});
        my $pos = $rating_hash{$score_sort[0]};
        my $text = '';
        if (defined $pos && $score_sort[0] != 0) {
            open my $fh, '>', \$text;
            draw_dimer($r->{f}, $r->{r}, $pos, $fh);
            close $fh;
        }
        $out = { dg => $dg + 0, diagram => $text };
    } elsif ($op eq 'orf') {
        my @orf = find_orf($r->{seq});
        $out = [ map { [ @$_ ] } @orf ];
    } elsif ($op eq 'cpg') {
        open my $null, '>', \my $discard;
        my $stdout = select $null;  # find_cpg prints debugging output
        my @cpg = find_cpg($r->{seq});
        select $stdout;
        $out = $cpg[0] == -1 ? [] : [ map { [ @$_ ] } @cpg ];
    } elsif ($op eq 'pcr') {
        $out = design($r, \&primer_window, 0);
    } elsif ($op eq 'bisulphite') {
        $out = design($r, \&bisul_window, 1);
    } elsif ($op eq 'qpcr') {
        %vars = %{ $r->{vars} };
        @intron_exon_bounds = @{ $r->{bounds} };
        my $seq = uc $r->{seq};
        $reverse = 0;
        @PF = primer_window($seq, length $seq);
        $pfkeys = @PF;
        $reverse = 1;
        @PR = primer_window(complement(scalar reverse $seq), length $seq);
        $qpcr_flag = 1;
        calc_amplicon(length $seq);
        $out = pairs_out();
    } elsif ($op eq 'sequencing') {
        %vars = %{ $r->{vars} };
        my $seq = uc $r->{seq};
        $reverse = 0;
        @PF = primer_window($seq, length $seq);
        $pfkeys = @PF;
        calc_seq_primers();
        $out = [ map { [ @$_ ] } @primer_pairs ];
    }
    print $json->encode({ result => $out }), "\n";
}
"""


def extract_subs(source: str) -> str:
    """Return the source of the subroutines listed in SUBS."""
    starts = [(m.start(), m.group(1)) for m in re.finditer(r"^sub (\w+)", source, re.M)]
    chunks = []
    for (start, name), (end, _) in zip(starts, [*starts[1:], (len(source), "")], strict=True):
        if name in SUBS:
            chunks.append(source[start:end])
    code = "\n".join(chunks)
    # "for $x qw(...)" is no longer valid syntax in modern Perl.
    return re.sub(r"for (\$\w+) qw\(([^)]*)\)", r"for \1 (qw(\2))", code)


def main() -> None:
    raw = Path(sys.argv[1]).read_bytes().decode("latin-1").replace("\r", "")
    sys.stdout.write(PRELUDE + extract_subs(raw) + DRIVER)


if __name__ == "__main__":
    main()
