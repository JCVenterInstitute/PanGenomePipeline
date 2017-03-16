#!/usr/bin/env perl
use strict;
use warnings;

use FileHandle;
use Getopt::Long;
use List::Util qw[min max];
#use Data::Dump 'dump'; 
my $properties_dir = "/usr/local/projdata/0695/projects/GENPROP/FTP/";
my $prop_def = "$properties_dir"."PROP_DEF";
my $prop_step = "$properties_dir"."PROP_STEP";
my $step_ev_link = "$properties_dir"."STEP_EV_LINK";
my $hmm_dir = "/usr/local/projdata/0695/projects/GENPROP/pfam_and_tigrfam/current/HMM/";
my $proj_name = "DEFAULT";
my $debug = 0;
my ( $seqs, $property, $list, $all, $eval_order );
GetOptions ( 'seqs=s' => \$seqs,
			'property=s' => \$property,
			'list=s' => \$list,
			'all' => \$all,
			'name=s' => \$proj_name,
			'eval_order=s' => \$eval_order,
            'debug=i' => \$debug,
        );
my %prop_def = ();
my %properties = ();
my %def_step = ();
my %steps = ();
my %step_ev = ();
my %ev = ();
my %results = ();
my @prop_list = ();
`curl ftp://ftp.jcvi.org/pub/data/TIGRFAMs/GEN_PROP/PROP_DEF_WITHOUT_DESCRIPTION_FIELD.TABLE > $prop_def` unless (-e $prop_def);
`curl ftp://ftp.jcvi.org/pub/data/TIGRFAMs/GEN_PROP/PROP_STEP.TABLE > $prop_step` unless (-e $prop_step);
`curl ftp://ftp.jcvi.org/pub/data/TIGRFAMs/GEN_PROP/STEP_EV_LINK.TABLE > $step_ev_link` unless (-e $step_ev_link);
open (DEF, "$prop_def");
open (RESULTS, ">>GP_results_$proj_name");
open (LONGFORM, ">>LONGFORM_REPORT_$proj_name");
open (TABLE, ">>TABLE_$proj_name");
open (SHORTFORM, ">>SUMMARY_FILE_$proj_name");

my $input = 0;
if ($all)
{
	$input++;
}
if ($list)
{
	$input++;
}
if ($property)
{
	$input++;
}
if ($input == 2)
{
	die "Only one form of input please\n";
}
if ($input > 2)
{
	die "Only one form of input please\n";
}
if (!$input)
{
	die "Please provide a form of input\n";
}

if ($list)
{
	open (LIST, "$list");
	while ( my $line = <LIST> )
	{
		chomp $line;
        my @input_list;
		push (@input_list, $line);
		@prop_list = @input_list;
	}
}

if ($eval_order && $all) 
{
	open (EVAL, "$eval_order");
	while ( my $line = <EVAL> )
	{
		my @split_line = split(/\s+/, $line);
		push (@prop_list, $split_line[0])
	}
	close EVAL;
}
if ($property)
{
	push (@prop_list, $property);
}
while ( my $line = <DEF> )
{
	my @split_line = split(/\t/, $line);															#splits up the prop def file lines by tab 				#in an array attached to the property name, and in a temp array for the line
	if ($split_line[2] eq "CATEGORY" || $split_line[2] eq "ROOT" || $split_line[2] eq "SUMMARY" ) {next;}
	$properties{$split_line[0]}[0] = $split_line[0];      								#store step_link????????
	$properties{$split_line[0]}[1] = $split_line[1]; 									#store name of genprop
	$properties{$split_line[0]}[2] = $split_line[2]; 									#store type of genprop
	$properties{$split_line[0]}[3] = $split_line[3]; 									#store GenProp acc
	$prop_def{$split_line[3]} = $split_line[0];
	if (!$eval_order && $all) {push(@prop_list, $split_line[3]);}	
}
close DEF;

open (STEP, "$prop_step");
while ( my $line = <STEP> )		#go through step table line by line
{

	my @split_line = split(/\t/, $line);															#split up the line of step table by tab
		if ($split_line[0] =~ /^\d/)		#if 2 number in each line of step table = step_link??????
		{
			$steps{$split_line[0]}[0] = $split_line[0];											#store 1st number									
			$steps{$split_line[0]}[1] = $split_line[1];											#store 2nd number
			$steps{$split_line[0]}[2] = $split_line[2];											#store 3rd number/word
			$steps{$split_line[0]}[3] = $split_line[3];											#store name
			$steps{$split_line[0]}[4] = $split_line[4];											#store next number
			push(@{$def_step{$split_line[1]}}, $split_line[0]);
		}
}
close STEP;	

open (EV, "$step_ev_link");
while ( my $line = <EV> )																			#go through it line by line
{
	my @split_line = split(/\t/, $line);															#split the line up by tabs
	#if ($split_line[1] == @{$steps{$identifier}[0]}[$k])										#if the 2nd number in the line equals the step identifier	
	$ev{$split_line[0]}[0] = $split_line[0];																				#store the 1st number
	$ev{$split_line[0]}[1] = $split_line[1];											#store 2nd number
	$ev{$split_line[0]}[2] = $split_line[2];											#store 3rd number/word
	$ev{$split_line[0]}[3] = $split_line[3];																			#store the type
	push(@{$step_ev{$split_line[1]}}, $split_line[0]);
}
close EV;

print TABLE ("acc\tdescription\tproperty_value\tstep_number\tstep_name\tstep_value\trequired?\tevidence_type\tevidence_name\tbest_HMM_hit\tHMM_hit_count\tHMM_command\n");

my $still_missing = 1;
while ($still_missing)
{
	$still_missing = check_results();
}

#####################################
sub check_results
{
	my $miss = 0;
	for my $acc (@prop_list)
	{
		#print "in check results: $acc\n";
		if (!defined($results{$acc}))
		{
			try_property($acc);
			if (!defined($results{$acc}))
			{
				if ($debug) {print "$acc is still missing\n";}
				$miss++;
				if ($debug) {print "value of missing: $miss\n";}
			}
			else
			{
				next;
			}
		}
		else
		{
			next;
		}
	}
	if ($debug) {print "MISSING is $miss <-----------------------------------\n";}
	return $miss;
}		
####################################

sub try_property
{
	my $acc = shift;
	my $acc_id = $prop_def{$acc};
	if ($debug) {print "prop is: $acc my acc_id is $acc_id\n";}
	my ( $answered, $result ) = evaluate_property($acc_id);      # <----- changed from 2 parameters to 1
	if ($answered)
	{
		$results{$acc} = $result;
		print_result($acc);
		if ($debug) {print "my result is $result\n";}
		return 1;
	}
	else
	{
		if ($debug) {print "cannot evaluate $result\n";}
		return 0;
	}
}		

###########################

sub evaluate_property
{
	my $prop_id = shift;
	my $found = 0;
	my $missing = 0;
    my ( $step_answered, $step_result );
	foreach (@{$def_step{$prop_id}})																#sort the hash keys by step number
	{
		my $identifier = $_;
		if ($debug) {print "step_id = $identifier\n";}
		if ($steps{$identifier}[4])													# check if required
			{
				my $required = 1;
				if ($debug) {print ("this step is required\n");}
			}
			else
			{
				my $required = 0;
				if ($debug) {print ("this step is not required\n");}
			}
			($step_answered, $step_result) = evaluate_step($identifier);
			if ($step_answered)
			{
				if ($step_result)
				{
					$steps{$identifier}[5] = "yes";
				}
				else
				{	
					$steps{$identifier}[5] = "no";
				}	
			}
			else
			{
				return (0,0);
			}
																		# We only need to report missing/found when a rule is required; missing/found are used for *property* evaluation
			if ($step_result)		#not required steps can turn no into partial								
			{
				$found = 1;
			}
			elsif ((!$step_result) && $steps{$identifier}[4]) #only required steps can turn yes into partial
			{
				$missing = 1;
			}
	}			
	#if ($debug) {print "missing:$missing\tfound:$found\n";}
	if (!$found)
	{
		return ($step_answered, 0);		# no required components were yes, thus the property is NOT FOUND (0)
	}
	elsif ($missing)
	{
		return ($step_answered, 1);		# at least one required component was a no (but at least one was found), thus property is PARTIAL (1) 
	}
	else
	{
		return ($step_answered, 2);		# no required components are missing (but at lesat one was found), thus property is YES (2)
	}
}

####################################################################################################

sub evaluate_step
{
	my $step_id = shift;
	my $succeed = 0;
	$some_evidence = 0;
	if ($debug) {print "I'm evaluating $step_id\n";}
	if (!defined($step_ev{$step_id}))
	{
		print ("\tthis step has no evidence\n");
	}
	else
	{
		$some_evidence = 1;
		foreach (@{$step_ev{$step_id}})
		{
			my $ev_id = $_;
			if (!$ev_id) {print RESULTS ("there is no ev_id\n"); next;}
			if ($debug) {print ("\tevidence is $ev_id\n");}
			if ($debug) {print "\tev type is: $ev{$ev_id}[3]\n";}
			if ($ev{$ev_id}[3] eq "RULE_BASE")																		#if the type is RULE_BASE
			{
				print RESULTS ("\tUNABLE TO EVALUATE EVIDENCE TYPE RULE_BASE. MARKING STEP AS NOT FOUND\n");					#warn that RULE_BASE cannot be evaluated at this time
				next;
			}	#end if RULE_BASE
			if (($ev{$ev_id}[3] eq "HMM") || ($ev{$ev_id}[3] eq "HMM-CLUST"))									#if the type is HMM-CLUST
			{
				print RESULTS ("\tSEARCHING HMM $ev{$ev_id}[2]\n");
				if ($ev{$ev_id}[3] eq "HMM-CLUST")
				{
					print RESULTS ("\tWARNING: DISTANCE CHECKING NOT CURRENTLY IMPLEMENTED. ONLY CHECKING HMM HIT\n");
				}
				$hmm_path = "$hmm_dir"."$ev{$ev_id}[2]".".HMM";
				if (! -e "$hmm_path")
				{
					print RESULTS ("WARNING: UNABLE TO FIND HMM $hmm_path\n MARKING STEP AS NOT FOUND\n");
					$ev{$ev_id}[4] = "BAD_HMM";														#store the first hit
					$ev{$ev_id}[5] = "BAD_HMM";				#store the number of hits
					next;
				}
				else
				{
					evaluate_hmm($hmm_path);
					#print"\tevaluating\n";
					$result = results_from_hmm("GP_temporary_HMM.tmp", $ev_id);
					if ($result)
					{
						$succeed++;
						#print "success1 = $succeed\n";
						next;
					}
					else 
					{
						next;
					}
				}	
			}	#end if HMM/HMM-CLUST
			if ($ev{$ev_id}[3] eq "GENPROP")																			#if the type is GenProp
			{
				if (!defined$prop_def{$ev{$ev_id}[2]})
				{
					print ("this step is dependent on non-existent genome property $ev{$ev_id}[2]. Marking as absent.\n");
					return (1,0);
				}
				if (defined$results{$ev{$ev_id}[2]})																	#if the results are known for the mentioned property
				{
					if ($results{$ev{$ev_id}[2]} == 0)					#NO is treated as NO
					{
						next;
					}
					elsif ($results{$ev{$ev_id}[2]} == 2  || ($results{$ev{$ev_id}[2]} == 1))			#YES and PARTIAL are treated as YES 
					{
						$succeed++;
						#print "success = $succeed\n";
						next;
					}
				}
				else
				{
					return (0,0);																							#unknown is unknown
				}
			}
		}
	}
	#print "SUCCESS = $succeed and evidence = $some_evidence\n";
	if ($succeed)
	{
		#print "step = success\n";
		return (1,1);
	}
	elsif ($some_evidence)
	{
		#print "evidence but false\n";
		return (1,0);
	}
	else
	{
		#print ("This step had no evidence; we are marking it as evaluated and false\n");
		return (1,0);
	}
}

#############################################################################################

sub evaluate_hmm
{
#	if (-e "GP_temporary_HMM.tmp") {`rm GP_temporary_HMM.tmp`;}
if (-s "GP_temporary_HMM.tmp" || -z "GP_temporary_HMM.tmp") {`rm GP_temporary_HMM.tmp`;}
	`hmm3search --noali --cut_tc --domtblout GP_temporary_HMM.tmp $hmm_path $seqs`;
}

##############################################################################################

sub results_from_hmm
{
	$source = shift;
	$ev_num = shift;
	$got_one = 0;
	open (HMMTABLE, "$source") || die "cannot read domtblout file\n";
	while ($line = <HMMTABLE>)
	{
		$first_hit = '';
		if (!($line =~ /^#/))																	#ignore the header line
		{
				$got_one=1;
				chomp $line;
		        @split_line = split(/\s+/, $line);
				$first_hit = $line; 
				last;                                       
		}
	}
	close HMMTABLE;
	if ($got_one) {$ev{$ev_num}[4] = $first_hit; $ev{$ev_num}[6] = "hmm3search --noali --cut_tc --domtblout GP_temporary_HMM.tmp $hmm_path $seqs";}														#store the first hit
	elsif (!$got_one) {$ev{$ev_num}[4] =  "<NO_HITS>";}
	$ev{$ev_num}[5] = `grep -Pv "^#" GP_temporary_HMM.tmp | wc -l`;				#store the number of hits
	chomp ($ev{$ev_num}[5]);
	if (!($ev{$ev_num}[5])) {$ev{$ev_num}[5] = 0;}
	return $got_one;
}

##############################################################################################

sub print_result
{
	my $prop = shift;
	my $prop_num = $prop_def{$prop};
	$printed_result = 0;
	$noev = 1;
	if ($debug) {print "PRINT: prop: $prop\n";}
	print LONGFORM ("PROPERTY: $prop\n$properties{$prop_num}[1]\n\n");
	print SHORTFORM ("$prop\t$properties{$prop_num}[1]\t");
	foreach (@{$def_step{$prop_num}})
	{
		$step_p = $_;
		if ($debug) {print "PRINT: step: $step_p\n";}
		print LONGFORM (".\tSTEP NUMBER: $steps{$step_p}[2]\n.\tSTEP NAME: $steps{$step_p}[3]\n");
		$already_printed = 0;
		foreach (@{$step_ev{$step_p}})
		{
			$noev = 0;
			$ev_p = $_;
			if ($debug) {print "PRINT: ev: $ev_p\n";}
			print TABLE ("$prop\t$properties{$prop_num}[1]\t");
			if ($results{$prop} == 2)
			{
				print TABLE ("YES\t");
				if (!$printed_result){print SHORTFORM ("YES\n"); $printed_result=1;}
			}
			if ($results{$prop} == 1)
			{
				print TABLE ("PARTIAL\t");
				if (!$printed_result){print SHORTFORM ("PARTIAL\n"); $printed_result=1;}
			}
			if ($results{$prop} == 0)
			{
				print TABLE ("NO\t");
				if (!$printed_result){print SHORTFORM ("NO\n"); $printed_result=1;}
			}
			print TABLE ("$steps{$step_p}[2]\t$steps{$step_p}[3]\t$steps{$step_p}[5]\t");
			if (@{$steps{$step_p}}[4])
			{
				print TABLE ("required\t");
				if (!$already_printed)
				{
					print LONGFORM (".\t.\trequired\n");
					$already_printed = 1;
				}
			}
			else
			{
				print TABLE ("not required\t");
				if (!$already_printed)
				{
					print LONGFORM (".\t.\tnot required\n");
					$already_printed = 1;
				}
			}
			print TABLE ("$ev{$ev_p}[3]\t$ev{$ev_p}[2]\t");
			if ($ev{$ev_p}[3] eq "HMM" || $ev{$ev_p}[3] eq "HMM-CLUST")
			{
				if ($ev{$ev_p}[5] ne "0")
				{
					@split_line = split(/\s+/, $ev{$ev_p}[4]);
					print TABLE ("$split_line[0]\t$ev{$ev_p}[5]\t$ev{$ev_p}[6]\n");
				}
				else
				{
					print TABLE ("$ev{$ev_p}[4]\t$ev{$ev_p}[5]\tN/A\n");
				}
				print LONGFORM (".\t.\t$ev{$ev_p}[3]: $ev{$ev_p}[2]\n.\t.\tBEST HIT: $ev{$ev_p}[4]\n.\t.\tHIT COUNT: $ev{$ev_p}[5]\n");
			}
			elsif ($ev{$ev_p}[3] eq "GENPROP")
			{
				print TABLE ("N/A\tN/A\tN/A\n");
				print LONGFORM (".\t.\t$ev{$ev_p}[3]: $ev{$ev_p}[2]\n");
			}		
			else 
			{
				print TABLE ("N/A\tN/A\tN/A\n");
				print LONGFORM (".\t.\t$ev{$ev_p}[3]\n");
			}
		}
		print LONGFORM (".\tSTEP RESULT: $steps{$step_p}[5]\n\n");
	}
	if ($noev)
	{
		print SHORTFORM ("NO\n");
	}
	if ($results{$prop} == 2)
	{
		print LONGFORM ("RESULT: YES\n\n//\n\n ");
	}
	if ($results{$prop} == 1)
	{
		print LONGFORM ("RESULT: PARTIAL\n\n//\n\n ");
	}
	if ($results{$prop} == 0)
	{
		print LONGFORM ("RESULT: NO\n\n//\n\n");
	}
}	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	



















