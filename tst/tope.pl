#!/usr/bin/perl

sub kldiv;

$ln2 = log(2);

open(F,"-|","./topetest",@ARGV);
while(<F>) {
	chomp;
	if (/Exact Result/) {
		while(<F>) {
			chomp;
			last if (/^\s*$/);
			my @p = split /\s+/, $_;
			push @exact, [ @p ];
		}
	} elsif (/Tope Result/) {
		my $i = 0;
		while(<F>) {
			chomp;
			last if (/^\s*$/);
			my @p = split /\s+/, $_;
			my $kl = kldiv($exact[$i],\@p);
			print "($kl) $_\n";
			$i++;
		}
	}
}


sub kldiv {
	my ($p1,$p2) = @_;
	my $ret = 0;
	for(my $i=0;$i<@$p1;$i++) {
		if ($p1->[$i] > 0) {
			$ret += $p1->[$i]*log($p1->[$i]/$p2->[$i])/$ln2;
		}
	}
	return $ret;
}
