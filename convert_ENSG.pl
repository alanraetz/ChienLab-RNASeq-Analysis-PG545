use strict;

my $input_file = $ARGV[0];

if ( ! $input_file ) { print "USAGE: $_ <DESeq2 CSV output file>\n"; exit; }

# my $ENSG_lookup = getEnsemblNames();
my ($ENSG_lookup,$name_lookup) = getHGNCNames();

open(IN,"$input_file") or die "Unable to open file $input_file\n"; 

my $prefix = $input_file; $prefix =~ s/\.csv|\.res\Z//;

my $named = $prefix . '_named.csv';

open(NAMED,">$named") or die "Unable to open -named file\n";
# process header
my($ENSG,$baseMean,$log2fold,$lfcES,$stat,$pvalue,$pAdj) = split(',',<IN>);
print NAMED "\"GeneName\",$baseMean,$log2fold,$lfcES,$stat,$pvalue,$pAdj";

my (%upreg,%downreg);
my $up_count = 0;
while ( my($ENSG,$baseMean,$log2fold,$lfcES,$stat,$pvalue,$pAdj) = split(',',<IN>) ) {
     
    $ENSG =~ s/\"//g;
    $ENSG =~ s/\.\d+\Z//; # remove suffix
    chomp $pAdj;
    if ( $pAdj eq 'NA' ) { next; }

    my $name = $ENSG_lookup->{$ENSG} || $ENSG;

    print NAMED "$name,$baseMean,$log2fold,$lfcES,$stat,$pvalue,$pAdj\n";

}

exit;


##############  END OF MAIN PROGRAM  ###################


sub getEnsemblNames {

    # ensembl.org => Biomart (select link at top of page) 
    #   => human genes => select only 2 attributes: 'Gene stable ID' and 'Gene name', click "Results"
    #   => select CSV output, click "Go", download file to local directory
    open(ENSG,"./ENSG-to-names.csv") or die "Unable to open Ensembl gene name lookup file ENSG-to-names.csv\n";
    my (%lookup,%name2ensg);
    while ( my($ensg,$name) = split(/,/,<ENSG>) ) {
        chomp($name);
        $ensg =~ s/\.\d+\Z//; # remove suffix
        $lookup{$ensg} = $name;
        $name2ensg{$name} = $ensg;
    }
    return (\%lookup,\%name2ensg); # reference to hash
}

#
#   https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_pub_ensembl_id&hgnc_dbtag=on&order_by=gd_pub_ensembl_id&format=text&submit=submit
#
sub getHGNCNames {

    open(HGNC,"./HGNC-to-names.txt") or die "Unable to open HGNC gene name lookup file ENSG-to-names.csv\n";
    my $header = <HGNC>; # first line
    my %lookup;
    while ( my($hgnc,$name,$ensg) = split(/\t/,<HGNC>) ) {
        chomp($ensg);
        $lookup{$ensg} = $name;
    }
    return \%lookup; # reference to hash

}
