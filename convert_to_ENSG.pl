use strict;

my $input_file = $ARGV[0];

if ( ! $input_file ) { print "USAGE: $_ <gene name list (one gene per row)>\n"; exit; }

my ($ENSG_lookup,$name_lookup) = getEnsemblNames();

my $hgnc_lookup = getHGNC_ENSG_id();

open(IN,"$input_file") or die "Unable to open file $input_file\n"; 

my $prefix = $input_file; $prefix =~ s/\.csv|\.txt\Z//;

my $named = $prefix . '_ENSG.csv';

open(NAMED,">$named") or die "Unable to open -named file\n";

my (%upreg,%downreg);
my $up_count = 0;

while ( my $geneName = <IN> ) {
     
    chomp $geneName;

    my $ensg = $name_lookup->{$geneName} || $hgnc_lookup->{$geneName};

    if ( $ensg eq '' ) { print "Could not find ENSG ID for $geneName!\n" ; $ensg = $geneName; }

    print NAMED "$ensg\n";

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
sub getHGNC_ENSG_id {

    open(HGNC,"./HGNC-to-names.txt") or die "Unable to open HGNC gene name lookup file ENSG-to-names.csv\n";
    my $header = <HGNC>; # first line
    my %lookup;
    while ( my($hgnc,$name,$ensg) = split(/\t/,<HGNC>) ) {
        chomp($ensg);
        $lookup{$name} = $ensg;
    }
    return \%lookup; # reference to hash

}
