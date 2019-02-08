#!/usr/bin/perl -w
use strict;
#This code was written by L. Sayavedra on Feb 2017

#change the line bellow to specify the path of a folder containing the subfolders with the Macsyfinder results
my @files = `ls ~/SecretionSystem/*/macsyfinder.out`;
my $out = "SecretionSystemSummaryTransporters.csv";
print @files;

my (%trans, $tr, $count, $genome, $name);
$genome =0; #counter for the genome being read
#print "Transporter \t$name\n";
foreach my $file (@files){
print "I am in genome $genome \n";
	open (MACSY, $file) or die;
	
	while (my $line = <MACSY>){
		if ($line =~ /.*\/(\S+)\.faa/){
			$name =$1;
			$trans {'Transporter'} [$genome]= $name;
		}
		
		if ($line =~ /(\S+)\t(\S+)/){
			$tr =$1;
			$count=$2;
			if (exists ($trans {$tr} [$genome])){  ###some systems, like VirB4 are detected with Conj, but not with other DB
				if ($count > $trans {$tr} [$genome]){
					$trans {$tr} [$genome]= $count;
				}
			}
		else {
				$trans {$tr} [$genome]= $count;
		}
		}
		
	}
	
	close (MACSY);
	$genome ++;
}

#print to a file the results
open (OUT, '>', $out) or die;
foreach my $amino (sort {$b cmp $a} (keys(%trans))) {#goes key by key
		 print "$amino,";
		 print OUT "$amino,"; 
			 foreach my $codon ( @{ $trans{$amino} } ) {      #returns the codons associated to a key 
					 print OUT "$codon,";
			} 
					 print OUT "\n"; 
}

close (OUT);
