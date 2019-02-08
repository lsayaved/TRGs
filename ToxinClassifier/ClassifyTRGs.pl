#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::DB::Fasta;

=pod

=head1 NAME

Classify potential TRG's


=head1 DESCRIPTION
This script will try to annotate the class to which the potential TRG belongs

=head1 AUTHOR

Lizbeth sayavedra


=head1 LENGUAGE

PERL

=head1 INPUT file
Needed: Output file of Conserved Protein Domain Family


=head1 OUTPUT Files

= file  with annotations

=head1 OPTIONS

-domains file with domains
-CPD output from CPD
-o, -output output file

=head1 HOW TO RUN THE PROGRAM

perl ClassifyTRGs.pl  -domains DomainsClassifier.txt -CPD ToxinSeacrh2_hitdata2.txt -out dummy.txt -faa fasta file with protein sequences
perl ClassifyTRGs.pl -domains DomainsClassifier.txt -CPD ToxinSearchV3_hitdata.txt -faa PossibleToxRelativesV2.fa  -out Test_out.txt
perl ClassifyTRGs.pl -domains DomainsClassifier.txt -CPD AllGenesCloseRelatives_hitdata.txt -faa ../SpeciesTree/fnnGenomes/Genes2Search.faa -out AllgenesTRGs -genids GenomeTags.txt
perl ClassifyTRGs.pl -domains DomainsClassifier.txt -CPD Genes2SearchPossibleTRGs_hitdata.txt -faa Genes2SearchPossibleTRGs.faa.fas  -out Test_out3 -genids GenomeTags.txt


=cut

###########Put in variables what the user is introducing#######
my %Argumentos = %{&ReadArguments};
my $file_domains =$Argumentos {'domains'};
my $CPD = $Argumentos{'CPD'};
my $output = $Argumentos{'out'};
my $seq_file = $Argumentos{'faa'};
#my $genids = $Argumentos {'genids'};
my $outputsum = $output ."_summ.txt";
my $outputdom = $output ."_dom.txt";
$output =$output. "pergenome.txt";
###############


###Main part

my $in; 
my $h=1;
my $TRG;
my $dom;
my $j=0; ##number of domains that have been classified
my %domains;
%domains = (domain => ['MARTX', 'RTX', 'RHS', 'RTXactivator']);
my %queries;
#%queries = (Queries => ['MARTX', 'RTX', 'RHS', 'RTXactivator', 'SUM', 'Class', 'Annotation']);
my %lengths;

print "Reading the domain scores (output from CD-Search Tool)... \n";
read_domain_scores (); 

open (CPD, $CPD) or die("Cannot open file: $CPD\n");
open (DOM, ">$outputdom");
print "Populating tables... \n";
while (<CPD>){
    if ($_ =~ /(\S+).*- >(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)/){
        my $query_ID =$1;
        my $Annotation = $2;
        $dom =$8; #this variable should have the short name of the domain that is in the main list ###full output is $8, short output from CPD is $9
        #print $dom . "\n";
	if (!exists ($queries{$query_ID})){  ##if the query already had at least one domain scored
            if (!exists ($domains{$dom})){
              push (@{$queries{$query_ID}},'0', '0', '0', '0', '0', 'X', $Annotation);
            print DOM "\n$query_ID\t$dom\t"; ##not sure about this
	    }
            else{
            print DOM "\n$query_ID\t$dom\t";
	    push (@{$queries{$1}}, $domains{$dom}[0], $domains{$dom}[1], $domains{$dom}[2], $domains{$dom}[3], '0', 'X', $Annotation);
            #print "$dom \t$queries{$query_ID}[0]\t$queries{$query_ID}[1]\t$queries{$query_ID}[2]\n";
            } 
        }
        else {
            if (!exists ($domains{$dom})){ ##if the domain doesn't exist don't add anything
		print DOM "$dom\t";  ## this will print to the domain file also the domains that I dont inlcude in my list 
            }    
            else {  ##if the domain exists, then the score increases for the deifned category
            print DOM "$dom\t";
            $queries{$query_ID}[0] = $queries{$query_ID}[0] + $domains{$dom}[0];
            $queries{$query_ID}[1] = $queries{$query_ID}[1] + $domains{$dom}[1];
            $queries{$query_ID}[2] = $queries{$query_ID}[2] + $domains{$dom}[2];
            ###test what domain is producing false positives in the RHS
            #print $dom. "\n";
	    $queries{$query_ID}[3] = $queries{$query_ID}[3] + $domains{$dom}[3];
            #print "$dom \t$queries{$query_ID}[0]\t$queries{$query_ID}[1]\t$queries{$query_ID}[2]\n";  
            }
        }
    }
}

close (CPD);
close (DOM);

print "Scoring the results... \n";

faa_length();
prescore ();
score();
print "Printing output in file $output \n";
print_output ();
#summary();

my $tag;
my $genome;
my %gentags;
my $dummy;
my $MARTX; my $RTX; my $RHS; my $fastaID; my $class; my $RTXactivator; my $Query;

###Subrutines
sub faa_length{
my $ID;
my $protein_length;
my $fasta_id;
my $seq;
$in  = Bio::SeqIO->new(-file => "$seq_file",
                       -format => 'Fasta');
    while ( $seq = $in->next_seq() ) {
        $ID ="Q#".$h;
        $fasta_id = $seq->id;
        $protein_length = $seq -> length; 
        push (@{$lengths{$ID}}, $fasta_id, $protein_length); ##first element of the array will contain identifier, second the length (1)
        $h++;
        #print "$ID\t $fasta_id\t$protein_length \n";
    }
}

sub prescore {
    foreach  $TRG (sort keys %lengths) {
        
        if (($lengths{$TRG}[1] > 1000) && ($queries{$TRG}[0]>=1) && ($queries{$TRG}[0]> $queries{$TRG}[2] )  ){ ##if protein is long, and the score of MARTX is bigger than 0, and the score of MARTX is lower than RHS
            if (exists ($queries{$TRG})){
                $queries{$TRG}[0] = $queries{$TRG}[0] + 100;
                print "$TRG\t$queries{$TRG}[0]\n";
            }
            else {
                print "Conflict\n\n";
            }
        }
        #elsif (($lengths{$TRG}[1] > 2000) && ($queries{$TRG}[0]>=1) && ($queries{$TRG}[0]> $queries{$TRG}[2] )   ){
        #    $queries{$TRG}[0] = $queries{$TRG}[0] + 30;
        #}
        #elsif (($lengths{$TRG}[1] < 1500) && ($queries{$TRG}[1]>=1) && ($queries{$TRG}[1]> $queries{$TRG}[2] )   ){
        #    $queries{$TRG}[1] = $queries{$TRG}[1] + 30;
        #}
    }
    
}

sub score {
##Which one had the best score?
#Sum of all scores along the rows, if the score is = 0 then there is no associated class
    foreach  $TRG (sort keys %queries) {
    $queries{$TRG}[4] = $queries{$TRG}[0] + $queries{$TRG}[1]  + $queries{$TRG}[2] ;
    }
    
    foreach $TRG (sort keys %queries) {
        if ($queries{$TRG}[4]==0){
            $queries{$TRG}[5] = "Unclassified";
          #  prescore();
        }
        
        
        elsif (($queries{$TRG}[0] == $queries{$TRG}[1]) && ($queries{$TRG}[0] > 0)  && ($queries{$TRG}[0] >$queries{$TRG}[2])){
            if (($lengths{$TRG}[1] < 2000 ) && ($queries{$TRG}[1]>=1) ){#&& ($queries{$TRG}[1]> $queries{$TRG}[2] )   ){
            $queries{$TRG}[5] = "RTX";
            print $TRG . "RTX";
	    }
            else {    
            $queries{$TRG}[5] = "RTX/MARTX";
            }
        }
        
        if (($queries{$TRG}[2] > $queries{$TRG}[1]) && ($queries{$TRG}[2] > $queries{$TRG}[0]) && ($queries{$TRG}[2] > 0.5) ){ ##added the >0.5 to avoid false positives
            $queries{$TRG}[5] = "RHS";
        }
        if (($queries{$TRG}[1] > $queries{$TRG}[0]) &&  ($queries{$TRG}[1]>=1) && ($queries{$TRG}[1] > $queries{$TRG}[2]) ){  # &&
            $queries{$TRG}[5] = "RTX";
	    
        }
        if (($queries{$TRG}[0] > $queries{$TRG}[1]) && ($queries{$TRG}[0] > $queries{$TRG}[2])  && ($queries{$TRG}[0]>=1)){  ##change the minimum value to 1 (it was with 0.07)
            $queries{$TRG}[5] = "MARTX";
        }
	if ($queries{$TRG}[3] > 0) {
            $queries{$TRG}[5] = "RTXactivator";
        }

    }

}



sub print_output {
#Print the result
open (RESULT, ">$output") or die("Cannot create file: $output\n");
print RESULT "Query\tMARTX\tRTX\tRHS\tRTXactivator\tSUM\tClass\tAnnotation\tType\tID\tLength";
foreach my $TRG (sort keys %lengths) {
    print RESULT "\n$TRG";
    if (exists ($queries{$TRG})){
        foreach (@{$queries{$TRG}}) {
            print RESULT "\t$_";
        }
        print RESULT "\t$lengths{$TRG}[0]\t$lengths{$TRG}[1]";
    }
    else {
        print RESULT "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$lengths{$TRG}[0]\t$lengths{$TRG}[1]"
        
    }
}

close (RESULT);
}

sub summary {
    #open (GENID, $genids);
    #while (<GENID>){
    #    if ($_ =~ /(\S+)\t(\S+)/){
    #        $genome=$1;
    #        $tag=$2;
    #        #$gentags{$tag}=$genome;  ##the key is the short name that is allways present in the identifiers of the faa sequence, the key is the genome name
    #        push (@{$gentags{$tag}}, $genome, '0', '0', '0', '0'); 
    #    }
    #}
    #close (GENID);
    
    open (READOUT, "$output") or die("Cannot open file: $output\n");
    open (OUT2, ">$outputsum") or die("Cannot open file: $output\n");
    print OUT2 "Query\tGenome\tGenID\tClass\tMARTX\tRTX\tRHS\tRTXactivator\n";
    while (<READOUT>){
	
        if ($_=~ /(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(.*)\t(.*)\t(\S+)/){ 
	    #print "$2\t$3\t$4\t$5\t$6\t$7\t$8\n";
	    $Query = $1;
	    $MARTX = $2;
            $RTX = $3;
            $RHS = $4;
	    $RTXactivator = $5;
            $class=$7;
            $fastaID = $9;  ###might be $10, depending on the format of the CPD file
            foreach my $key (keys %gentags) {
                if ($fastaID =~ /$key/){  ##format is str, substr, str contains substr
                #print "I work with the short name $key that should be contained in $fastaID\n";
                print OUT2 "$Query\t$gentags{$key}[0]\t$fastaID\t$class\t$MARTX\t$RTX\t$RHS\t$RTXactivator\n";   
                }
            }
        }
    }
            
    
    close (READOUT);
    close (OUT2);

}


sub read_domain_scores{
####put in an array all the identifiers that are read in a file
open (FILE, $file_domains);
while (my $LINE = <FILE>){
    if ($LINE =~ /(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(.*)/){
     push (@{$domains{$1}}, $2, $3, $4, $5);
    }
}
close (FILE);
}




sub ReadArguments {
    my %opts = ();
    GetOptions (\%opts,
		'domains|d=s',
                'CPD=s',
		'out=s',
                'faa=s',
                'genids=s',
		'help|h');
    
    if(($opts{'help'}) || (scalar(keys(%opts))==0)){  #se llama a la subrutina ayuda en caso de que se escriba help, o en caso de que el usuario no haya metido ningun parametro
		&PrintHelp()  
	}
	
	return \%opts;    #me devuelve una referencia a la tabla opts
}
	

####Subrutina que despliega la ayuda en caso de que el usuario no haya escrito las opciones necesarias, se llama desde el pod
sub PrintHelp {  
    system "pod2text -c $0";
    exit()
}

