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
        &qsubMemThreads
        &qsubJava
        &qsubJavaMaster
        &jobNative
   );
}

#EXCEPTIONS
#my %excepts=();
#$excepts{"MEM_THREADS"}{"REALIGNMENT"}=1;
#$excepts{"MEM_THREADS"}{"CALLING"}=1;

sub generic(@){
    my ($opt,$function)=@_;
    my $qsub = "qsub -P ".$$opt{CLUSTER_PROJECT}." -pe threaded ".$$opt{$function."_THREADS"}." -q ".$$opt{$function."_QUEUE"};
    return ($qsub);
}

sub qsubTemplate(@){
    my ($opt,$function)=@_;
    my $qsub = &generic($opt,$function)." -m a -M ".$$opt{MAIL}." -R ".$$opt{CLUSTER_RESERVATION}." -l h_rt=".$$opt{$function."_TIME"}.",h_vmem=".$$opt{$function."_MEM"}."G";
    return ($qsub)
}
sub qsubMemThreads(@){
    my ($opt,$function)=@_;
    my $h_vmem= ($$opt{$function."_MEM"} * $$opt{$function."_THREADS"})."G";
    my $qsub = &generic($opt,$function)." -m a -M ".$$opt{MAIL}." -R ".$$opt{CLUSTER_RESERVATION}." -l h_rt=".$$opt{$function."_TIME"}.",h_vmem=".$h_vmem;
    return ($qsub)
}

sub qsubJavaMaster(@){
    my ($opt,$function)=@_;
    my $h_vmem = (4 + $$opt{$function."_MASTERTHREADS"} * $$opt{$function."_MEM"})."G";
    my $qsub = &generic($opt,$function)." -m a -M ".$$opt{MAIL}." -R ".$$opt{CLUSTER_RESERVATION}." -l h_rt=".$$opt{$function."_TIME"}.",h_vmem=".$h_vmem;
    return ($qsub)
}

sub qsubJava(@){
    my ($opt,$function)=@_;
    my $h_vmem = (4 + $$opt{$function."_THREADS"} * $$opt{$function."_MEM"})."G";
    my $qsub = &generic($opt,$function)." -m a -M ".$$opt{MAIL}." -R ".$$opt{CLUSTER_RESERVATION}." -l h_rt=".$$opt{$function."_TIME"}.",h_vmem=".$h_vmem;
    return ($qsub)
}


sub jobNative(@){
    my ($opt,$function)=@_;
    my $h_vmem = (4 + $$opt{$function."_MASTERTHREADS"} * $$opt{$function."_MEM"})."G";
    my $qsub = &generic($opt,$function)." -l h_rt=".$$opt{$function."_TIME"}.",h_vmem=".$h_vmem;
    return ($qsub)
}

1;
