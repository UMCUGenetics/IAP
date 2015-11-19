#
# Copyright (c) Hindrik H.D. Kerstens and
# UMC Utrecht. All rights reserved
#


package illumina_sge;
use strict;
use warnings;

BEGIN {

    use Exporter;

    our @ISA = ('Exporter');

    our @EXPORT = qw (
        &qsubTemplate
   );
}




sub qsubTemplate(@){
    my ($opt,$function)=@_;
    my $qsub = "qsub -P ".$$opt{CLUSTER_PROJECT}. " -m a -M ".$$opt{MAIL}." -pe threaded ".$$opt{$function."_THREADS"}." -R ".$$opt{CLUSTER_RESERVATION}.
        " -l h_rt=".$$opt{$function."_TIME"}.",h_vmem=".$$opt{$function."_MEM"}."G";
    return ($qsub)
}



1;
