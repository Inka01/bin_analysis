#!/usr/bin/env perl
###############################################################################
#
#    Bin_summary.pl
#
#    Reads in scaffold data and summarizes into table
#
#    Copyright (C) Inka Vanwonterghem
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

#core Perl modules
use Getopt::Long;
use Carp;
#use Math;

#CPAN modules
use Bio::SeqIO;
use Data::Dumper;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
printAtStart();
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################

# read files

#format:
#   contig_ID   cov1    cov2    cov3    cov4    cov5    ... cov[i]
#   contig_521  0.220   1.445   0       0.445   5.500   ... 18.941
my $coverages = openRead($global_options->{'coverages'});

#format:
#   BIN_ID  completeness    contamination
#   BIN_1   90.09           2.7
#   BIN_50  54.59           0.9
my $check = openRead($global_options->{'checkm'});

#folder containing all the BIN_*.fa files
my $dir = $global_options->{'bin_folder'};
opendir (DIR, $dir) or die "can't open the directory";
my @bins = readdir DIR;
close DIR;

# write to a file
my $out = openWrite($global_options->{'summary'});

#global variables
my %cov = ();
my %check = ();
my @column = ();
my $max_index = -1;

#read the coverage file and create a hash of arrays to contain all the coverage information per contig
foreach my $contigline (<$coverages>) {
    next if ($contigline =~ /^#/);
    next if ($contigline =~ /^$/);
    chomp $contigline;
    my @info = split(/\t/,$contigline);
    my @tmp_array = ();
    if($max_index == -1) {
        $max_index = $#info-1;
    }
    foreach my $j (1..$max_index+1) {
        push(@tmp_array, $info[$j]);
    }
    $cov{$info[0]} = \@tmp_array;
}

#read the checkM file and create a hash of arrays to contain all the completeness/contamination information per bin
foreach my $binline (<$check>) {
    next if ($binline =~ /^#/);
    next if ($binline =~ /^$/);
    chomp $binline;
    my @info = split(/\t/,$binline);
    my @tmparray = ($info[1],$info[2]);
    $check{$info[0]} = \@tmparray;
}

foreach my $i (1..$max_index+1) {
    my $string = sprintf("avecov%d",$i);
    push (@column,$string);
}
print $out "BIN_ID\tcontigs\tlength\tN50\tcompleteness\tcontamination\t".join("\t",@column)."\n";

#go through each of the BIN*.fa files, use bioperl to obtain the seq_id and seq_length per contig
foreach my $bin (@bins) {
    next if ($bin !~ /.fa$/);
    my @cov_array = ();
    my $N50 = 0;
    my $avecov;
    my @avecov_array;
    my $bin_ID;
    my $completeness;
    my $contamination;
    $bin_ID = $bin;
    $bin_ID =~ s/.fa//;
    #create an array to save all the lengths
    my @lengths = ();
    #create a scalar to store the total length
    my $lengthsum = 0;
    #create an array to store the coverage in order to calculate the average for the total bin
    foreach my $i (0..$max_index) {
        push (@cov_array,0);
    }
    #create a scalar to store the number of contigs in the bin
    my $contigcount = 0;
    #create scalars to store the completeness and contamination
    $completeness = 0.00;
    $contamination = 0.00;
    #obtain the seq_id (contig name) and seq_length for each sequence in the fasta file, add the length to the array, obtain the coverage profiles for this contig, and add to contigcount
    my $seqio = Bio::SeqIO->new(-file => $bin,'-format' => 'Fasta');
    while (my $seq = $seqio->next_seq) {
        my $seq_string = $seq->seq;
        my $seq_id = $seq->id;
        #print "$seq_id";
        my $seq_len = length($seq_string);
        #my @tmp_array = \@lengths;
        push (@lengths, $seq_len);
        $lengthsum += $seq_len;
        $contigcount += 1;
        $seq_id =~ s/_contig/contig/;
        if (exists $cov{$seq_id}) {
            foreach my $i (0..$max_index) {
                $cov_array[$i] += $cov{$seq_id}->[$i];
            }
        }
    }
    #print Dumper(@cov_array);
    #now we need to calculate the N50 value using the array of lengths
    my $cutoff = $lengthsum/2;
    my @sortedlengths = sort { $b <=> $a } @lengths;
    my $sum = 0;
    my $tempN50 = $sortedlengths[0];
    $N50 = 0;
    #my $i = 0;
    #while ($i <= $#sortedlengths) {
        #$sum += $i;
        #last if ($sum >= $cutoff);
        #$i += 1;
    #}
    #$N50 = $sortedlengths[$i];
    foreach (@sortedlengths) {
        my $value = $sum + $_;
        if ($value > $cutoff) {
            $N50 = $tempN50;
            last;
        }
        $sum += $_;
        $tempN50 = $_;
    }
    #obtain the completeness and contamination from the checkM file
        if (exists $check{$bin_ID}) {
            $completeness = $check{$bin_ID}->[0];
            $contamination = $check{$bin_ID}->[1];
        }
    #print $completeness;
    #print $contamination;
    #now we need to print all the information to a file
    $avecov = 0.00;
    foreach my $j (0..$max_index) {
        $avecov = $cov_array[$j]/$contigcount;
        push (@avecov_array,$avecov);
    }

    print $out "$bin\t$contigcount\t$lengthsum\t$N50\t$completeness\t$contamination\t".join("\t",@avecov_array)."\n";
}

close ($out);
close ($coverages);
close ($check);
exit;


######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    # "long|single<type>"
    #
    # :i number (integer)
    # :f number (decimal)
    # :s string
    # +  flag
    #
    my @standard_options = ( "help|h+", "coverages|c:s", "checkm|m:s", "bin_folder|b:s", "summary|o:s");
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    #if(!exists $options{''} ) { printParamError (""); }
    if(!exists $options{'coverages'} ) { printParamError ("You need to supply the coverage file"); }
    if(!exists $options{'checkm'} ) { printParamError ("You need to supply the checkm results file"); }
    if(!exists $options{'bin_folder'} ) { printParamError ("You need to supply the folder containing the *.fa bin files"); }
    if(!exists $options{'summary'} ) { printParamError ("You need to supply an output file name"); }
    
    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #
    my ($error) = @_;
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub createDefault
{
    #-----
    # Set default values for parameters
    #
    my ($option_name, $default_value) = @_;
    if(not exists $global_options->{$option_name})
    {
      $global_options->{$option_name} = $default_value;
    }
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          },
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;

    isCommandInPath($cmd, $failure_type);

    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};

    my $cmd_str = $cmd . " " . $param_str;

    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to
    # checkAndRunCommand
    #
    my $ref = shift;

    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}


######################################################################
# MISC

sub printAtStart {
print<<"EOF";
----------------------------------------------------------------
 $0
 Copyright (C) Inka Vanwonterghem

 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
----------------------------------------------------------------
EOF
}

__DATA__

=head1 NAME

    Bin_summary.pl

=head1 COPYRIGHT

   copyright (C) Inka Vanwonterghem

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

    This script takes the sspace output (after parsing with script Jason -
    Scaffolding_binned_contigs_Jason) and creates two new files: one summary
    of all the inter-intra-unbinned links for each bin (total and percentage),
    one summary of all the links for each bin

=head1 SYNOPSIS

    Bin_summary.pl -c <coverages.txt> -m <checkm.txt> -b <BIN_folder> -o <summary.txt> [-help|h]

    -coverages -c               input file containing the coverage for each contig in each sample
    -checkm -m                  input file containing the checkm results per bin (contamination and completeness)
    -bin_folder -b              folder containing all the BIN*.fa files (make sure you are in this folder)
    -summary -o                 output file containing the summary for each bin
    [-help -h]                  Displays basic usage information

=cut


