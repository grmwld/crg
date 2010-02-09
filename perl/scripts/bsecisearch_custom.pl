#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Std;

# Command line options
my %opts;
getopts('i:s:f:a:o:', \%opts);

# Predefined values
my @alphabet = ('A','C','G','U');
my $original_file = 'UGAend_sequence';
my $input_file = 'RNAfold.result';
my $output_file = 'bsecis_containing_sequences';
my $orf_file = 'orf_sequences';
my $rna_file = 'rna_sequences';
my $gc_file = 'gc_content.out';

my $energy_limit = -7.5;
my $min_pos_apical_loop = 19;
my $max_pos_apical_loop = 50;
my $min_length_apical_loop = 3;
my $max_length_apical_loop = 9;
my $min_length_upper_stem = 3;
my $max_length_upper_stem = 16;

my $min_window_length = 39;
my $max_window_length = 80;

my $min_length_before_secis = 2;
my $min_length_after_secis = 0;
my $min_length_orf = 15;


my $genome_file; 
my $annotation;
my $secis_sequence;
my $orf;
my $orf_length;
my $protein;

my $uga_count = 0;
my @contents = ( );
my %genome = ( );
my %sequences = ( );
my %secondary_structures = ( );
my $free_energy;

my $secis_start_position = 0;
my $secis_end_position = 0;
my $genome_id;
my @genome_seq = ('','');

my $genome_length;

my $char;
my $secis_length;
my $sequence_id;
my $sequence;
my $ss;

my $extended_length = 3;

my $num = 0;

my $output_dir = '.';
if ($opts{'o'}) {
    $output_dir = $opts{o}
}

open(TEMP, ">$output_dir/$output_file") or die "Can't open file $output_file.\n";
open(ORF, ">$output_dir/$orf_file") or die "Can't open file $orf_file.\n";
open(RNA, ">$output_dir/$rna_file") or die "Can't open file $rna_file.\n";
open (GC_CON, ">$output_dir/$gc_file") or die "Can't open file.\n";

# command line arguments parsing
my $dostrands;
my $dostrands_start;
my $doframes;
my $step;
# input file
if ($opts{'i'}) {
    $genome_file = $opts{'i'}
}
else
{
    print "Please input query genome/contig file name:\n";
    $genome_file = <STDIN>;
}

chomp $genome_file;

# strands to check. Default is both
if ($opts{'s'}) {
    $dostrands = int($opts{'s'}) - 1;
    $dostrands_start = $dostrands
}
else {
    $dostrands = 1;
    $dostrands_start = 0
}

# frames to check. Default is all (3)
if ($opts{'f'}) {
    $doframes = int($opts{'f'}) - 1;
    $step = 3
}
else {
    $doframes = 0;
    $step = 1;
}

print "Loading data......";
sleep (1);
@contents = get_file_data($genome_file);
%genome = extract_sequence_from_fasta_data(@contents);
print "OK!\n";
sleep (1);

foreach $genome_id (keys %genome) {
     
      $genome_seq[0] .= $genome{$genome_id};  
}

$genome_seq[0] = uc $genome_seq[0];

$genome_seq[0] =~ s/T/U/g;

$genome_seq[0] =~ s/[^ACGU]//ig;

$genome_length = length $genome_seq[0];

$genome_seq[1] = reverse $genome_seq[0];

$genome_seq[1] =~ tr/ACGU/UGCA/;


my @number_of_N;
my @frequency_of_N;

my $strand_num;
my $c;
my $genome_seq;

print "\n\n";
print "**************   Running bSECISFinding    ******************\n";
print "Genome length: ",$genome_length, " letters\n";

for (my $strand_num = $dostrands_start; $strand_num <= $dostrands; $strand_num++) {

      	$genome_seq = $genome_seq[$strand_num];

	@number_of_N = (0,0,0,0);
	@frequency_of_N = (0,0,0,0);

	for (my $pos = 0; $pos < $genome_length; $pos++) {
      
      		$c = substr($genome_seq, $pos,1);

      		if ( $c eq 'A') {
      		   $number_of_N[0]++;
      		} elsif ($c eq 'C') {
		   $number_of_N[1]++;
      		} elsif ($c eq 'G') {
		   $number_of_N[2]++;
      		} else {
        	   $number_of_N[3]++;
      		}
	}
	$genome_seq = '';

	for (my $i = 0; $i < 4; $i++) {
      		$frequency_of_N[$i] = $number_of_N[$i]/$genome_length;
	}

	if ($strand_num == 0) {
		print "The nucleotide frequencies in plus(+) strand are:\n";
		print GC_CON "The nucleotide frequencies in plus(+) strand are:\n";
	} else { 
		print "The nucleotide frequencies in minus(-) strand are:\n";
		print GC_CON "The nucleotide frequencies in minus(-) strand are:\n";
	}
		
	for (my $i = 0; $i < 4; $i++) {
      		print GC_CON $alphabet[$i]," -> ";
		printf GC_CON "%3.3f", $frequency_of_N[$i];
		print GC_CON "\t";
	}

	print GC_CON "\n";
}


$secis_start_position = 0;

print "\n\n";

for (my $strand_num = $dostrands_start; $strand_num <= $dostrands; $strand_num++) {

	$genome_seq = $genome_seq[$strand_num];

	if ($strand_num == 0) {
		print "\n**********************  Now searching the plus(+) strand.......\n\n";
	} else {
		print "\n**********************  Now searching the minus(-) strand.......\n\n";
	}

	for ($secis_start_position = $doframes; $secis_start_position < $genome_length - $min_window_length; $secis_start_position+=$step) {

    	 	$orf = '';
     	 	$secis_sequence = '';
     	 	$secis_length = 0;
   
     	 	if (($char = substr($genome_seq, $secis_start_position, 3)) ne 'UGA') {
		
                	next;
      		} else {
	
			$uga_count++;
			$secis_sequence .= $char;
		
			$secis_end_position = $secis_start_position + 3;
                	$secis_length = 3;

			### ORF constraints

			my ($orf_start_position, $orf_end_position, $upstream_stop_position);

			$upstream_stop_position = $secis_start_position;

			my $codon;
			my $aug_number = 0;

			while (($upstream_stop_position - 3 >= 0) && !(($codon = substr($genome_seq, $upstream_stop_position - 3, 3)) =~ /UAA|UAG|UGA/gm)) 
			{
				
				$upstream_stop_position = $upstream_stop_position - 3;

				$orf = $codon.$orf;

				if (($codon eq 'AUG') || ($codon eq 'GUG')) {
					$aug_number = 1;
					last;
				}

			}
			
			if ($upstream_stop_position < 3) {
				next;
			}
			
			
			if ($aug_number == 0) {
				next;
			}

			if ((length $orf) < $min_length_before_secis) {
				next;
			}

			$orf_start_position = $upstream_stop_position;
			
	
	     	LOOP: if ($secis_end_position >= $genome_length)  
			{
				next;
			}

			if (($codon = substr($genome_seq, $secis_end_position, 3)) =~ /UAG|UAA|UGA/gm)
			{
				#next;
			}

			#while (!(($char = substr($genome_seq, $secis_end_position, 3)) =~ /UAG|UAA|UGA/gm)) 
			while ($char = substr($genome_seq, $secis_end_position, 3))
			{
				$secis_sequence .= $char;
                        	$secis_length += 3;
				$secis_end_position += 3;

                        	if ($secis_length >= $min_window_length) {
                        		last;
                        	}
                	}

			if (($secis_length < $min_window_length) || ($secis_length >= $max_window_length)) {
                  	      next;
                	}

			if ($strand_num == 0) {

				$annotation = '>SECIS element (Location: (+) '.($secis_start_position + 1).' - '.$secis_end_position.', SECIS Length = '.($secis_end_position -$secis_start_position).' nt)';
			} else {
				$annotation = '>SECIS element (Location: (-) '.($genome_length - $secis_start_position).' - '.($genome_length - $secis_end_position + 1).', SECIS Length = '.($secis_end_position -$secis_start_position).' nt)';
			}

			@contents = ( );

            		open(ORIGIN, ">$original_file") or die "Can't open file $original_file.\n";
			print ORIGIN $annotation;
			print ORIGIN "\n";
			print ORIGIN $secis_sequence;
			
			close ORIGIN;

			if (system("RNAfold <$original_file >$input_file") != 0) {
 				print "Can't predict secondary structure of original data";
      				#exit;
			}

			@contents = get_file_data($input_file);
			%sequences = extract_sequence_and_structure(@contents);
			$free_energy = get_free_energy(@contents);

			if ($free_energy >= $energy_limit) {
				goto LOOP;
			}
	
			@contents = ( );

			$sequence_id = (keys %sequences)[0];
			
			my @ss_segments = ( );
      			my $pos;
      			my ($pos1, $pos2);
      			my $length;
      			my $len1 = 0;
      			my $len2 = 0;
      
      			my ($apical_loop, $left_upper_stem, $right_upper_stem, $upstream_sequence, $downstream_sequence);
      			my ($apical_loop_start, $apical_loop_end);
      			my ($left_upper_stem_start, $left_upper_stem_end, $right_upper_stem_start, $right_upper_stem_end);
      			my ($upstream_sequence_end, $downstream_sequence_start);
      			
      			$sequence = $sequences{$sequence_id}->[0];
      			$ss = $sequences{$sequence_id}->[1];
      			$length = length $sequence;

      			if ($ss =~ /\d+.\d+/) {
          			next;
      			}
			
      			####    Step 1: searching for apical loop
			
			$pos = $min_pos_apical_loop;
			
			if ((($char = substr($ss, $pos, 1)) eq ')') || (($char = substr($ss, $pos - 1, 1)) eq ')') || (($char = substr($ss, $pos + 1, 1)) eq ')')) {
				while (($char = substr($ss, $pos, 1)) ne '(') {
					$pos++;
					
					if ($pos >= ((length $ss) - 1)) {
						goto LOOP;
					}
				}
				
			}

      			while (($char = substr($ss, $pos, 1)) ne ')') {
				$pos++;

				if ($pos >= ((length $ss) - 1)) {
					goto LOOP;
				}
      			}

			if ($pos >= ((length $ss) - 1)) {
				goto LOOP;
			}

      			$apical_loop_end = $pos-1;

      			while (($char = substr($ss, $pos, 1)) ne '(') {
           			$pos--;
      			}

      			$apical_loop_start = $pos+1;

				###  Apical loop length and composition constraints

			if (($apical_loop_start < $min_pos_apical_loop) || ($apical_loop_end > $max_pos_apical_loop)) {
				goto LOOP;
			}

			
      			$apical_loop = substr($sequence, $apical_loop_start, ($apical_loop_end - $apical_loop_start+1));


			my $apical_length = length $apical_loop;

			my $found_con1 = 0;
			my $found_con2 = 0;
			
			if (($apical_length == 3) && ($apical_loop =~ /GU|GG|GA/)) {
				$found_con1 = 1;
				
				
			} elsif ($apical_loop =~ /GU|GG|GA/) {
				if ((length $` <= 2) && (length $` <= ($apical_length/2 - 1))) {
					$found_con1 = 1;
				}

			} else {
				;
			}

			#if (($found_con1 == 0) && ($apical_loop =~ /GAA|AAA/))  {
				
			#	if ((length $` <= 1) && (length $` <= ($apical_length/2 - 1))) {
			#		$found_con2 = 1;
			#	}
			#}


			if ($found_con1 == 1) {

				if ((length $apical_loop < $min_length_apical_loop) ||(length $apical_loop > $max_length_apical_loop)) {
					goto LOOP;
				}
			} elsif ($found_con2 == 1) {
				
				if ((length $apical_loop < 3) || (length $apical_loop > 14)) {
					goto LOOP;
				}
			} else {
				goto LOOP;
			}
			
    			

      			####    Step 2: searching for upper stem

      			$left_upper_stem_end = $apical_loop_start - 1;
      			$right_upper_stem_start = $apical_loop_end + 1;
      			$pos1 = $left_upper_stem_end;
      			$pos2 = $right_upper_stem_start;

      			while (($char = substr($ss, $pos1, 1)) eq '(') {
	  			$pos1--;
           			$len1++;
      			}

      			while (($char = substr($ss, $pos2, 1)) eq ')') {
	   			$pos2++;
           			$len2++;
      			}

      			while ($len1 < $len2) {           
        			LINE1: while (($char = substr($ss, $pos1, 1)) eq '.') {
               
                			$pos1--;
           			}  
 
           			while (($char = substr($ss, $pos1, 1)) eq '(') {
                
                			$pos1--;
                			$len1++;
           			}

				if ((($char = substr($ss, $pos1, 1)) eq ')') && ($len1 != $len2)) {
					goto LOOP;
				}

           			if ($len1 == $len2) {
                			last;
           			} elsif ($len1 < $len2) {
                			goto LINE1;
					
           			} else {
               				 LINE2: while (($char = substr($ss, $pos2, 1)) eq '.') {
                      				$pos2++;
                			 }

                			while (($char = substr($ss, $pos2, 1)) eq ')') {
                      				$pos2++;
                      				$len2++;
                			}

					if ((($char = substr($ss, $pos2, 1)) eq '(') && ($len1 != $len2)) {
						goto LOOP;
					}

                			if ($len1 == $len2) {
                   				last;
                			} elsif ($len2 < $len1) {
                      				goto LINE2;
                			} else {
                      				goto LINE1;
                			}
          			}
      			}

      			while ($len1 > $len2) {
        			LINE3: while (($char = substr($ss, $pos2, 1)) eq '.') {
                				$pos2++;
           			}
           
           			while (($char = substr($ss, $pos2, 1)) eq ')') {
                			$pos2++;
                			$len2++;
           			}

				if ((($char = substr($ss, $pos2, 1)) eq '(') && ($len1 != $len2)) {
						goto LOOP;
					}

           			if ($len1 == $len2) {
                			last;
           			} elsif ($len1 > $len2) {
               	 			goto LINE3;
           			} else {
               			LINE4: while (($char = substr($ss, $pos1, 1)) eq '.') {
                      				$pos1--;
                			}
	
                			while (($char = substr($ss, $pos1, 1)) eq '(') {
                      				$pos1--;
                      				$len1++;
                			}

					if ((($char = substr($ss, $pos1, 1)) eq ')') && ($len1 != $len2)) {
						goto LOOP;
					}

                			if ($len1 == $len2) {
                   				last;
                			} elsif ($len1 < $len2) {
                      				goto LINE4;
                			} else {
                      				goto LINE3;
                			}
          			}
      			}

			if ($found_con1 == 1) {

				if (($len1 < $min_length_upper_stem) || ($len1 > $max_length_upper_stem) || ($len2 < $min_length_upper_stem) || ($len2 > $max_length_upper_stem)) {
				goto LOOP;
				}
			}
			
			

			if ($found_con2 == 1) {
				if (($len1 < 4) || ($len1 > 10) || ($len2 < 4) || ($len2 > 10)) {
				goto LOOP;
				}
			}

      			$left_upper_stem_start = $pos1 + 1;
      			$right_upper_stem_end = $pos2 - 1;
      			$left_upper_stem = substr($sequence, $left_upper_stem_start, ($left_upper_stem_end - $left_upper_stem_start + 1));
      			$right_upper_stem = substr($sequence, $right_upper_stem_start, ($right_upper_stem_end - $right_upper_stem_start + 1));

			if (($left_upper_stem eq '') || ($apical_loop eq '') || ($right_upper_stem eq '')) {
				
				goto LOOP;
			}

			
			####    Step 3: searching for upstream and downstream sequence
			
			$upstream_sequence_end = $left_upper_stem_start - 1;
			$downstream_sequence_start = $right_upper_stem_end + 1;


			$upstream_sequence = substr($sequence, 0, ($upstream_sequence_end - 0 + 1));

			my $s = substr($ss, 0, ($upstream_sequence_end - 0 + 1));
 
			#if ($s =~ /\)/) {goto LOOP;}

			$downstream_sequence = substr($sequence, $downstream_sequence_start, ((length($sequence) - 1) - $downstream_sequence_start + 1));

 
			$orf = $orf.' '.$secis_sequence.' ';
			$orf_end_position = $secis_end_position;

			while (($orf_end_position <= ($genome_length - 3)) && !(($codon = substr($genome_seq, $orf_end_position, 3)) =~ /UAA|UAG|UGA/gm)) {
				$orf_end_position = $orf_end_position + 3;
				$orf = $orf.$codon;
			}

			
			if ($orf_end_position > ($genome_length - 3)) {
				#$orf_end_position = $genome_length - 1;
				next;
			}


			if (($orf_end_position - $secis_end_position) < $min_length_after_secis) {
				next;
			}

			if (($codon = substr($genome_seq, $orf_end_position, 3)) =~ /UAA|UAG|UGA/gm) {
				$orf_end_position = $orf_end_position + 3;
				$orf = $orf.$codon;
			}else {
				next;
			}


			$orf_length = (length $orf) - 2;

			if ($orf_length < $min_length_orf) {
				next;
			}


			$num++;

			#chop $sequence_id;
			#$sequence_id .= ', Free energy = '.$free_energy.' kcal/mol)';
			print "No. ",$num,"\t";
			print $sequence_id, "\n";

      			print TEMP $sequence_id, "\n";
			print TEMP $upstream_sequence, "\t\t\t";
      			print TEMP $left_upper_stem, "\t";
      			print TEMP $apical_loop, "\t";
      			print TEMP $right_upper_stem, "\t\t\t";
			print TEMP $downstream_sequence;
      			print TEMP "\n";

			print TEMP substr($ss, 0, ($upstream_sequence_end - 0 + 1)), "\t\t\t";
			print TEMP substr($ss, $left_upper_stem_start, ($left_upper_stem_end - $left_upper_stem_start + 1)), "\t";
			print TEMP substr($ss, $apical_loop_start, ($apical_loop_end - $apical_loop_start + 1)), "\t";
			print TEMP substr($ss, $right_upper_stem_start, ($right_upper_stem_end - $right_upper_stem_start + 1)), "\t\t\t";
			print TEMP substr($ss, $downstream_sequence_start, ((length($sequence) - 1) - $downstream_sequence_start + 1));
			
	
			print TEMP "\n";

			print ORF $sequence_id, ", ORF length = ", $orf_length, " nt (", ($orf_length/3), "aa)\n";
			print RNA $sequence_id, ", ORF length = ", $orf_length, " nt (", ($orf_length/3), "aa)\n";
			print RNA $orf,"\n";

			$orf =~ s/\s//gm;
			$protein = translate_frame($orf, 1);
			chop($protein);
			
			print ORF $protein,"\n";

		}
	}
}


close TEMP;
close ORF;

print "\nThere are ", $uga_count, " UGA codons in the query genome (two strands).\n";
print "There are ", $num, " possible SECIS elements in the query genome.\n\n";

sleep (1);

print "Running bSECISProfile algorithm............";

system("rm $input_file");
system("rm $original_file");
system("rm *.ps");

print "OK!!\n";

#### Subroutines #####

sub get_file_data {

    my ($filename) = @_;
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}

sub extract_sequence_from_fasta_data {

    my (@fasta_file_data) = @_;
    my $annotation = '';
    my $sequence = '';
   
    my %anno_seqinfo = ();
 
    foreach my $line (@fasta_file_data) {

        if ($line =~ /^\s*$/) {
            next;

        } elsif($line =~ /^\s*#/) {
            next;
    
        } elsif($line =~ /^>/) {
            $sequence='';
	    $annotation = $line;
            $annotation =~ s/[\t\n\r\f]//g;
            next;

        } else {
	    $sequence .= $line;
            $sequence =~ s/\s//g;
	    $sequence =~ s/\d//g;
            $anno_seqinfo{$annotation}=$sequence;
	} 
    }

    return %anno_seqinfo;
}


sub extract_sequence_and_structure {
    my (@fasta_file_data) = @_;
    my $annotation = '';
    my $sequence = '';
    my $secondary_structure = '';
    my $energy;
   
    my %anno_seqinfo = ();
 
    foreach my $line (@fasta_file_data) {

        if ($line =~ /^\s*$/) {
            next;

        } elsif($line =~ /^\s*#/) {
            next;
    
        } elsif($line =~ /^>/) {
            $sequence='';

	    $annotation = $line;
            $annotation =~ s/[\t\n\r\f]//g;
            next;

        } elsif($line =~ /^[AGCU]/) {
	    $sequence .= $line;
            $sequence =~ s/\s//g;
            $anno_seqinfo{$annotation}->[0]=$sequence;
            
	} else {
            $secondary_structure = $line;

	    if ( $secondary_structure =~ /\s+\( *(-\d+.\d+)\)/gm ) {
		$energy = $1;
	    }

	    $secondary_structure =~ s/\s+\( *-\d+.\d+\)//g;
            
            $anno_seqinfo{$annotation}->[1]=$secondary_structure;
	  
        }
    }

    return %anno_seqinfo;
}


sub get_free_energy {

    my (@data) = @_;
    my $energy;

    foreach my $line (@data) {

        if ($line =~ /\s+\( *(-{0,1}\d+.\d+)\)/gm ) {
	    $energy = $1;
	} else {
	    next;
	}
    }

    return $energy;
}

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'UCA' => 'S',    # Serine
    'UCC' => 'S',    # Serine
    'UCG' => 'S',    # Serine
    'UCU' => 'S',    # Serine
    'UUC' => 'F',    # Phenylalanine
    'UUU' => 'F',    # Phenylalanine
    'UUA' => 'L',    # Leucine
    'UUG' => 'L',    # Leucine
    'UAC' => 'Y',    # Tyrosine
    'UAU' => 'Y',    # Tyrosine
    'UAA' => 'x',    # Stop
    'UAG' => 'x',    # Stop
    'UGC' => 'C',    # Cysteine
    'UGU' => 'C',    # Cysteine
    'UGA' => 'u',    # Sec
    'UGG' => 'W',    # Tryptophan
    'CUA' => 'L',    # Leucine
    'CUC' => 'L',    # Leucine
    'CUG' => 'L',    # Leucine
    'CUU' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCU' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAU' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGU' => 'R',    # Arginine
    'AUA' => 'I',    # Isoleucine
    'AUC' => 'I',    # Isoleucine
    'AUU' => 'I',    # Isoleucine
    'AUG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACU' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAU' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGU' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GUA' => 'V',    # Valine
    'GUC' => 'V',    # Valine
    'GUG' => 'V',    # Valine
    'GUU' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCU' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAU' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGU' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
            return '*';
    }
}

sub rna2peptide {

    my($rna) = @_;

    use strict;
    use warnings;
   

    my $protein = '';

    for(my $i=0; $i < (length($rna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($rna,$i,3) );
    }

    return $protein;
}

sub translate_frame {

    my($seq, $start, $end) = @_;

    my $protein;

    unless($end) {
        $end = length($seq);
    }

        return rna2peptide ( substr ( $seq, $start - 1, $end -$start + 1) );
}
