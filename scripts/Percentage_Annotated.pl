#! /usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use List::Util qw(sum);
use List::MoreUtils qw(any);
use Data::Dumper;


#all global variables needed
my $gb_file = $ARGV[0];
my $protein_id_CDS = "";
my $fasta_len = 0;
my $start = 0;
my $end = 0;

#open .gb file
my $seqio_object = Bio::SeqIO->new(-file => $gb_file, -format => "genbank");
my $entry = $seqio_object->next_seq;
my $id = $entry->id;




for my $feat_object ($entry->get_SeqFeatures) {
        if ($feat_object->primary_tag eq "CDS"){
		#loop for the first CDS: extracting name, locus tag and coordinates
		if ($protein_id_CDS eq ""){
			$protein_id_CDS = ($feat_object->get_tag_values("protein_id"))[0];
			my @locations = $feat_object->location->each_Location;
			my $first = $locations[0];
			$start = $first -> start();
			my $last = $locations[-1];
			$end = $last -> end();
		}
		
		#loop for all next CDSs
		else{
			#write
         		my $path = "CDSs.txt";
         		open(my $f, ">>", $path);
         		my $str = $id. "\t". $protein_id_CDS. "\t". $start. "\t". $end;
         		say $f $str;
         		close $f;

 			#set names for the current CDS
			$protein_id_CDS = ($feat_object->get_tag_values("protein_id"))[0];
			my @locations = $feat_object->location->each_Location;
			my $first = $locations[0];
			$start = $first -> start();
			my $last = $locations[-1];
			$end = $last -> end();
			
		}


                
        }
           

#if mat_peptide is found under the CDS (and has the proper locus tag), add its length to "annotated"
        if ($feat_object->primary_tag eq "mat_peptide"){
			my @locations = $feat_object->location->each_Location;
			my $first = $locations[0];
			my $start = $first -> start();
			my $last = $locations[-1];
			my $end = $last -> end();
			my $name = ($feat_object->get_tag_values("product"))[0];
			#write
         		my $path = "mat_peptides.txt";
         		open(my $f, ">>", $path);
         		my $str = $id. "\t". $protein_id_CDS. "\t". $name. "\t". $start. "\t". $end;
         		say $f $str;
         		close $f;

        }

}



