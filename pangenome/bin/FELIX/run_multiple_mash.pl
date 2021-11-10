#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd;
use File::Basename;

my $dirname = dirname(__FILE__);

my ($second_diff, $mash_file, $mash_ver, $name_add, $input, $out, $help, $cut_off, $run_phy, $tmp_path, $show_cnt, $also_run, $run_all_v_all, $is_read, $is_contig, $mat_only, $verb, $self, $kmer, $size, $trim_file, $is_reset, $version, $quiet, $input_file, $input_col, $id_col, $prev, $type_strain, $notable);
$tmp_path = getcwd;
my $MAX = 1000;
$input_col = 1;
$name_add = ".fasta";
$size = 1000;
GetOptions("second|w=s"=>\$second_diff, "add|a=s"=>\$name_add, "num|e=s"=>\$MAX, "phy|u=s"=>\$run_phy, "mash_file|m=s"=>\$mash_file, "mash_exec|M=s"=>\$mash_ver, "read|g"=>\$is_read, "all|j=s"=>\$run_all_v_all, "size|z=s"=>\$size, "trim|r=s"=>\$trim_file, "only_matrix|y"=>\$mat_only, "cnt|b"=>\$show_cnt, "kmer|k=s"=>\$kmer, "self|s"=>\$self, "notable|x"=>\$notable, "version|v=s"=>\$version, "quiet|q"=>\$quiet, "input|i=s"=>\$input, "out|o=s"=>\$out, "help|h|?"=>\$help, "cutoff|c=s"=>\$cut_off, "tmp_dir|t=s"=>\$tmp_path, "contig|t" =>\$is_contig, "reset|n" =>\$is_reset, "file|f=s"=>\$input_file, "col|l=s"=>\$input_col, "id|d=s"=>\$id_col, "prev|p=s"=>\$prev, "type_strain|T=s"=>\$type_strain);
if ($help || !(($mash_file || $run_all_v_all) && ($input || $input_file)))
{
	print STDERR "Runs MASH of a directory or fasta file against the MASH sketch\n\n";
	print STDERR "--------------------USAGE--------------------------------------\n";
	print STDERR "	-i input string containing directory with the fasta files, the fasta file, or an ls style search string to use (this has to be in \"quotes\")\n";
	print STDERR "	-f input file containing a list of the files to use. Default column to use is 0.\n";
	print STDERR "	-l column in the input file to use using 1-count. Default column to use is 1.\n";
	print STDERR "	-d column in the input to use as a name. Otherwise use the file id.\n";
	print STDERR "	-p column in the input to use as a previous identification. Not used otherwise\n";
	print STDERR "	-o output file prefix\n";
	print STDERR "	-c comma seperated cutoff values. Best hits within this value are combined. Default is 97,94.\n";
	print STDERR "	-j flag to run the input string/files as an all versus all. Ignores m. Should be the number of kmers in the sketch file\n";
	print STDERR "	-m mash scatch file to compare to\n";
	print STDERR "	-M mash executable\n";
	print STDERR "	-s remove self vs self output in the text file\n";
	print STDERR "	-t run files as a contigs\n";
	print STDERR "	-T type strain identifier to use for flitering and first row and column\n";
	print STDERR "	-v version of mash to use (1.1 or 2.0). Default of 2.0\n";
	print STDERR "	-q quiet mode\n";
	print STDERR "	-k kmer sketch size. Default is sketch file\n";	
	print STDERR "	-r trim file\n";	
	print STDERR "	-b show kmer counts in report table\n";
	print STDERR "	-n reset 0 to minimal\n";
	print STDERR "	-x does not make the table\n";
	print STDERR "	-y only makes the matrix\n";
	print STDERR "	-g files are reads. Incompatable with -t\n";
	print STDERR "	-e number of files to use per sub-run. Default is 1000.\n";
	print STDERR "	-a string to add to name id for mapping file. Default is '.fasta'\n";
	print STDERR "	-u string of the type to run phylogeny \n";
	print STDERR "	-w 2nd best hit ANI difference which to combine. Default is 0\n";
	print STDERR "	-? help\n";
		
	exit(0);
}

my @cut_off;
if ($cut_off)
{
	@cut_off = split ",", $cut_off;
}
else
{
	@cut_off = (97, 94);
	$cut_off = "97,94";
}
push @cut_off, 0;

my $add;
if ($size && !$run_all_v_all)
{
	$add .= "-s $size ";
}
if ($kmer && !$run_all_v_all)
{
	$add .= "-k $kmer ";
}


if (!$quiet)
{
	$verb = 1;
}

if (!$out)
{
	$out = $tmp_path . "/temp.";
}

my $name_file = $mash_file; $name_file =~ s/.msh/.txt/;
if ($run_all_v_all)
{
	$mat_only = 1;
}
if (!$run_all_v_all)
{
	if (!-e $mash_file && !-e $name_file )
	{
		print STDERR "Cannot locate $name_file and $mash_file. Please retry..\n";
		exit(0);
	}
	if ( !-e $name_file && !$mat_only)
	{
		print STDERR "Cannot locate name file... only making matrix..\n";
		$mat_only = 1;
	}
}
my %name;

if (!$mash_ver && ($version))
{
	if ($version ne "2.0" && $version ne "1.1")
	{
		print STDERR "MASH version $version not recognized. Using version 2.0 ...\n\n";
		$mash_ver = "$dirname/mash2";
	}
	else
	{
		if ($version eq "1.1")
		{
			$mash_ver = "$dirname/mash.sh";
		}
		else
		{
			$mash_ver = "$dirname/mash2";
		}
	}
}
if (!$mash_ver)
{
    $mash_ver = "$dirname/mash";
}
if (!$run_all_v_all && -e $name_file)
{
	open(my $fh, "<", $name_file);
	while(my $in = <$fh>)
	{
		while($in =~ /([^\n\r]+)/g)
		{	
			my @in = split "\t", $1;
			$name{$in[0] . $name_add}->{path} = $in[2];
			$name{$in[0] . $name_add}->{name} = $in[1];
		}		
	}
	close($fh);
}
my $arr;
my $refs;
my $dir_list;

my $list;
my $list;
my $out1;
my $dir;

my %trim_list;
if (-e $trim_file)
{
	open(my $t_file, "<", $trim_file);
	while(my $v = <$t_file>)
	{
		chomp($v); $v =~ /\A(\S+)/; $trim_list{$1} = 1; 
	}
}
if ($input)
{
	if ($input =~ /,/)
	{
		$input .= ",";
		$input =~ s/,/\n/g;
		$dir = $input;
	}
	else
	{
		if (-d $input)
		{
			$input .= "/*";
		}
		$dir = `ls $input`;
	}
}
my %new_names;
my %prev_val;
if ($input_file)
{
	if (-e $input_file)
	{
		open( my $fh, "<", $input_file); 
		while( my $v = <$fh> )
		{
			chomp($v); my @n = split "\t", $v; 
			if (-e $n[$input_col-1])
			{
				$dir .= $n[$input_col-1] . "\n";
			}
			if ($id_col)
			{
				$n[$input_col-1] =~ /([^\/\\]+)\Z/;
				$new_names{$1} = $n[$id_col-1];
			}
			if ($prev)
			{
				$n[$input_col-1] =~ /([^\/\\]+)\Z/;
				$prev_val{$1} = $n[$prev-1];
			}
		}
	}
}
my $cnt2;
my @out_list;
my $n = 0;
my $q_tax;
my $cnts;
my $temp_mash_id;
if ($dir)
{
	if ($run_all_v_all)
	{
		my $dir1 = "";
		while ($dir =~ /([^\n\r]+)([\n\r])/g)
		{
			my $tmp = $1; my $t1;
			$tmp =~ /([^\/]+)\Z/; $t1 = $1; 
			if ($trim_file)
			{
				my $id1 = $1; 
				$id1 =~ /(\S+)\.([^\.]+)\Z/;
				my $id = $1; 
				if ($trim_list{$id})
				{
					$dir1 .= " $tmp";
				}
			}
			else
			{
				if (-e $tmp && !-d $tmp)
				{
					$dir1 .= " $tmp";
				}				
			}
		}
		my $k = 21;
		if ($kmer) { $k = $kmer; }
		$temp_mash_id = rand(100000);
		my $addk = "-k $k";
		my $add1 = " -s $run_all_v_all -o $temp_mash_id ";
		if ($is_contig && !$is_read)
		{
			$add1 .= "-i ";
		}
		if ($is_read)
		{
			$add1 .= "-r ";
		}
		my $hide = `$mash_ver sketch $addk $add1 $dir1 2>&1`;
		while ($hide =~ /WARNING/) 	
		{ 
			$hide =~ /a(\s)k-mer(\s)size(\s)of(\s)at(\s)least(\s)(\d+)/;
			print STDERR "kmer warning give for $kmer size. K of $7 is recommended. Re-running sketch with new kmer size....\n";
			$k = $7;
			$addk = "-k $k";
			$hide = `$mash_ver sketch $addk $add1 $dir1 2>&1`;
		}
		$add1 = $addk . $add1;
		if (! -e ($temp_mash_id . ".msh"))
		{
			warn("Cannot make mash sketch file... quitting...\n");
			exit(1);
		}
		$mash_file = $temp_mash_id . ".msh"
	}
	while ($dir =~ /([^\n\r]+)([\n\r])/g)
	{
		my $tmp = $1; my $t1;
		$cnt2++;
		if ($cnt2 == $MAX) { $cnt2 = 0; $n++; }
		$tmp =~ /([^\/]+)\Z/; $t1 = $1; 
		if ($trim_file)
		{
			my $id1 = $1; 
			$id1 =~ /(\S+)\.([^\.]+)\Z/;
				my $id = $1; 
				if ($trim_list{$id})
				{
					$out_list[$n] .= " $tmp";
				}
		}
		else
		{
			if (-e $tmp && !-d $tmp)
			{
				$out_list[$n] .= " $tmp";
			}				
		}
		if ($tmp =~ /\.msh\Z/)
		{
			my $qtxt = $tmp; $qtxt =~ s/msh\Z/txt/;
			if (-e $qtxt)
			{
				open(my $qfo, "<", $qtxt);
				while( my $v1 = <$qfo>)
				{
					chomp($v1); my @n_1 = split "\t", $v1;
					$q_tax->{$n_1[0] . $name_add} = clean_name($n_1[1],100);
					
				}
				close($qfo);
			}
		}
		$dir_list->{$t1} = $tmp;
	}
		
		for (my $i = 0; $i < scalar(@out_list); $i++)
		{
				$out1 = $out_list[$i];
				print STDERR "run $i... \n";
				if ($is_contig)
				{
					$list = `$mash_ver dist -i $add $mash_file $out1`;
				}
				else
				{
					if ($is_read)
					{
						$list = `$mash_ver dist -r $add $mash_file $out1`;
					
					}
					else
					{
						$list = `$mash_ver dist $add $mash_file $out1`;
					}
			
				}
		
				while ($list =~ /([^\n\r]+)([\n\r])/g)
				{
					my $in = $1; my @in = split "\t", $in;
					$in[0] =~ /([^\/]+)\Z/; my $ref = $1;
					$in[1] =~ /([^\/]+)\Z/; my $go = $1;
					if ($name{$ref} || $mat_only)
					{
						$arr->{$go}->{$ref} = sprintf("%.2f", 100*(1-$in[2]));
						$in[4] =~ /(\d+)\/(\d+)/; 
						$cnts->{$go}->{$ref} = $in[4];
						$refs->{$ref} = 1;
					}
				}
		}
	}
	else
	{
		print STDERR "Unable to find anything... quitting\n\n";
		exit(0);
	}

my @refs = keys(%$refs); #add a sort here? how are we sure this is the same sort as the row order
my $out_str = "ID";
my $out_str2;

my %next_list;
my $best_id;
my $cut_off_sp;# $cut_off_sp= 94;
my $oth_list;
my $oth_list_cnt;

my $oth_list2;
my $min_val;
my %cut_off_sp_val;

foreach my $a (@refs) { my $a1 = $a; if ($new_names{$a}) { $a1 = $new_names{$a};  } $out_str .= "\t$a1"; }  $out_str .= "\n";
$out_str2 = $out_str;
print STDERR scalar(keys(%$arr)), "\n";
my %best2_list;
foreach my $b (keys(%$arr))
{
	 my $a1 = $b; if ($new_names{$b}) { $a1 = $new_names{$b};  } 
	 if (!$notable)
	 {
		$out_str .= $a1; $out_str2 .= $a1; 
		$cut_off_sp_val{$b} = "NA";
	}
	 foreach my $a (@refs) 
	 {
	    if (!$notable)
		{
			$out_str .= "\t". $arr->{$b}->{$a};
			$out_str2 .= "\t". $cnts->{$b}->{$a};
			
	    }
		if ($arr->{$b}->{$a} > 0 && ( !$min_val || $min_val > $arr->{$b}->{$a}))
	    {
			$min_val = $arr->{$b}->{$a};
	    }
		for (my $i = 0; $i < scalar(@cut_off); $i++)
		{
			if ($arr->{$b}->{$a} > $cut_off[$i]) 
			{
				$next_list{$b}->[$i] .= "$a(". $cnts->{$b}->{$a}."),";
			}
		}
		if ($a ne $b || !$self)
		{
			my $hit;
			if ($second_diff >0)
			{
				if (!$best2_list{$b} || $best2_list{$b}->[0]->{id} < $arr->{$b}->{$a})
				{
					if ($best2_list{$b} && clean_name($name{$best2_list{$b}->[0]->{ref}}->{name}, $best2_list{$b}->[0]->{id}) ne clean_name($name{$a}->{name}, $arr->{$b}->{$a})) 
					{
						$best2_list{$b}->[1]->{id} = $best2_list{$b}->[0]->{id};
						$best2_list{$b}->[1]->{ref} = $best2_list{$b}->[0]->{ref};
					}
					$best2_list{$b}->[0]->{id} = $arr->{$b}->{$a};
					$best2_list{$b}->[0]->{ref} = $a;
				}
				else
				{
					if ($best2_list{$b}->[1]->{id} < $arr->{$b}->{$a} && clean_name($name{$best2_list{$b}->[0]->{ref}}->{name}, $best2_list{$b}->[0]->{id}) ne clean_name($name{$a}->{name}, $arr->{$b}->{$a}))
					{
						$best2_list{$b}->[1]->{id} = $arr->{$b}->{$a};
						$best2_list{$b}->[1]->{ref} = $a;
					}
				}
			}
			for (my $i = 0; $i < scalar(@cut_off); $i++)
			{
				if (!$hit)
				{
					if($arr->{$b}->{$a} > $cut_off[$i] &&  ($arr->{$b}->{$a} > $best_id->{$b}->{id} || ( $arr->{$b}->{$a} == $best_id->{$b}->{id} && length($name{$a}->{name}) > length($name{$best_id->{$b}}->{name}))))
					{
						$hit = 1;
						$cut_off_sp_val{$b} = $cut_off[$i];
						if ($name{$best_id->{$b}->{name}}->{name} ne $name{$a}->{name})
						{
							$best_id->{$b}->{name2} = $best_id->{$b}->{name}; 
							$best_id->{$b}->{best2} = $best_id->{$b}->{id};
							$best_id->{$b}->{cnt2} = $best_id->{$b}->{cnt};
						}	
			
						$best_id->{$b}->{name} = $a; 
						$best_id->{$b}->{id} = $arr->{$b}->{$a};
						$best_id->{$b}->{cnt} = $cnts->{$b}->{$a};
					}
					else
					{
						if($arr->{$b}->{$a} > $cut_off[$i] &&  $arr->{$b}->{$a} > $best_id->{$b}->{best2} && $name{$a}->{name} ne $name{$best_id->{$b}->{name}}->{name})
						{
							$hit = 1;
							$best_id->{$b}->{best2} = $arr->{$b}->{$a}; $best_id->{$b}->{name2} = $a; $best_id->{$b}->{cnt2} = $cnts->{$b}->{$a};
						}
					}
				}
			}
			my $cl = clean_name($name{$a}->{name}, 100); 
			if (!$oth_list->{$b}->{$cl} || $oth_list->{$b}->{$cl} < $arr->{$b}->{$a})
			{
				$oth_list->{$b}->{$cl} = $arr->{$b}->{$a};
				$oth_list_cnt->{$b}->{$cl} = $cnts->{$b}->{$a};
			}
			
		}
	 }
	if (!$notable)
	{
		
		$out_str .= "\n";
		$out_str2 .= "\n";
		
	}
}
if ($is_reset)
{
    print STDERR "Resertting all 0 ANI to non-zero minimal value: $min_val\n";
    my $z_cnt;
    $out_str = "";
	
    foreach my $a (@refs) { my $a1 = $a; if ($new_names{$a}) { $a1 = $new_names{$a};  } $out_str .= "\t$a1"; }  $out_str .= "\n";
    print STDERR scalar(keys(%$arr)), "\n";
    foreach my $b (keys(%$arr))
    {
	my $a1 = $b; if ($new_names{$b}) { $a1 = $new_names{$b}; } $out_str .= $a1;
	foreach my $a (@refs)
	{
	    if ($arr->{$b}->{$a} > 0)
	    {
			if (!$notable)
			{		
				$out_str .= "\t". $arr->{$b}->{$a};
			}
	    }
	    else
	    {
			$z_cnt++;
			if (!$notable)
			{		
			
			$out_str .= "\t" . $min_val;
			}
		}
	}
	if (!$notable)
	{		
		$out_str .= "\n";
	}
    }
    print STDERR "Found $z_cnt cells with 0 value...\n";
}

if (!$notable)
{
	if ($out)
	{
		open(my $fo, ">", $out . "mat"); print $fo $out_str; close($fo);
		open(my $fo, ">", $out . "cnt"); print $fo $out_str2; close($fo);
	
	}
	else
	{
		print $out_str;
	}
}

my $fo;
if (!$mat_only)
{
	if ($out)
	{
		open($fo, ">", $out . "txt"); 
	}
	else
	{
		$fo = *STDOUT;
	}
	my $Q1 = "Query_name";
	if ($q_tax) { $Q1 .= "\tQuery_Taxonomy"; }
	if ($also_run)
	{
		my $gani_prog_path = "/usr/local/devel/ANNOTATION/tclarke/bin/ANIcalculator";
	
		my $gani2 = "perl /home/tclarke/ANI.pl -bl /usr/local/packages/ncbi-blast/bin/blastall -fd /usr/local/packages/ncbi-blast/bin/formatdb -od $tmp_path/tmp/";
		print $fo "$Q1\tBest MASH Reference Name\tBest MASH Full Reference Organism Name\tBest MASH Clean Reference Organism Name\tBest MASH ANI\tBest JSpecies Reference\tBest JSpecies Full Reference Organism Name\tBest JSpecies Reference Organism Clean Name\tBest JSpecies Reference ANI\n";

		foreach my $new (keys(%next_list))
		{
			if ($verb)
			{
				my $qnew = $new; if ($q_tax) { if ($q_tax->{$new}) { $qnew .= "\t". $q_tax->{$new}; } else {$qnew .= "\tNA";}}
				print STDERR "$qnew\t", $best_id->{$new}->{name}, "\t", $name{$best_id->{$new}->{name}}->{name}, "\t", clean_name($name{$best_id->{$new}->{name}}->{name}, $best_id->{$new}->{id}), "\t", $best_id->{$new}->{id},"\t",$name{$best_id->{$new}->{name2}}->{name},"\t", $best_id->{$new}->{best2},"\n";
	
			}
			my $full_1 = $dir_list->{$new};
			my @good_refs = split ",", $next_list{$new};
			my $gani_best_id;

			foreach my $file (@good_refs)
			{
				my $full_2 = $name{$file}->{path};
				my $sys = "rm $tmp_path/ANIcalculator.*";
				#`$sys`;
				$sys = "rm -r $tmp_path/ani.blast.dir/";
				#`$sys`;
		
				#my $in = `$gani_prog_path -stdout -genome1fna $full_1 	-genome2fna $full_2 -outdir $tmp_path`;
				#$in =~ /([^\r\n]+)\Z/;
				#my $val = read_log_file($tmp_path . "/ANIcalculator.log", 1);
				my $in2 = `$gani2 -qr $full_1 -sb $full_2`;
				chomp($in2);
				my @val = split "\t", $in2;	
				print STDERR $in2, "\t",$name{$val[1]}->{name},"\n";
			
				if (!$gani_best_id || $val[2] > $gani_best_id->{id})
				{
					$gani_best_id->{id} = $val[2];
					$gani_best_id->{name} = $file;
				}
				$sys = "rm -r $tmp_path/tmp/";
				`$sys`;
		
			}
			my $nameID = $new;
			if ($new_names{$nameID}) { $nameID = $new_names{$nameID}; }
			my $qnew = $nameID; if ($q_tax) { if ($q_tax->{$new}) { $qnew .= "\t". $q_tax->{$new}; } else {$qnew .= "\tNA";}}
		
			print $fo "$qnew\t", $best_id->{$new}->{name}, "\t", $name{$best_id->{$new}->{name}}->{name}, "\t", clean_name($name{$best_id->{$new}->{name}}->{name}, $best_id->{$new}->{id}), "\t", $best_id->{$new}->{id},"\t",$gani_best_id->{name}, "\t", $name{$gani_best_id->{name}}->{name}, "\t", clean_name($name{$gani_best_id->{name}}->{name}, $gani_best_id->{id}), "\t", $gani_best_id->{id}, "\n";

		}
	}
	else
	{
		if (!%cut_off_sp_val)
		{
			print $fo "$Q1\tReference Name\tFull Reference Organism Name\tClean Reference Organism Name\tReference ANI\t2nd Best Organism\t2nd Best ANI\n";
			foreach my $a (keys(%$best_id))
			{
				my $out_id = $a;
				#if ($name{$a}->{name}) {$out_id = $name{$a}->{name}}
				if (!$best_id->{$a}->{id})
				{
					$best_id->{$a}->{id} = "NA";
				}
				my $nameID = $a;
				if ($new_names{$nameID}) { $nameID = $new_names{$nameID}; }
	
				my $qnew = $nameID; if ($q_tax) { if ($q_tax->{$a}) { $qnew .= "\t". $q_tax->{$a}; } else {$qnew .= "\tNA";}}

				print $fo "$qnew\t", $best_id->{$a}->{name}, "\t", $name{$best_id->{$a}->{name}}->{name}, "\t", clean_name($name{$best_id->{$a}->{name}}->{name}, $best_id->{$a}->{id}, $cut_off), "\t(", $best_id->{$a}->{id},",", $best_id->{$a}->{cnt},")\t",$name{$best_id->{$a}->{name2}}->{name},"\t", $best_id->{$a}->{best2},"\n";
			}
		}
		else
		{
			if ($prev)
			{
				print $fo "$Q1\tReference Name\tPrevious Name\tFull Reference Organism Name\tClean Reference Organism Name\tReference ANI\tCutoff Value\t# of Other Hits\tOther Hits\t2nd Best Organism\t2nd Best ANI\n";
			}
			else
			{
				print $fo "$Q1\tReference Name\tFull Reference Organism Name\tClean Reference Organism Name\tReference ANI\tCutoff Value\t# of Other Hits\tOther Hits\t2nd Best Organism\t2nd Best ANI\n";
			}
			foreach my $a1 (keys(%$best_id))
			{
				my $out_id = $a1;
				my @s = sort {$oth_list->{$a1}->{$b} <=> $oth_list->{$a1}->{$a}} keys(%{$oth_list->{$a1}});
				#if ($name{$a}->{name}) {$out_id = $name{$a}->{name}}
				my $c = 1; my $tmp; 
				if ($cut_off_sp_val{$a1} > 0) 
				{ 
					while($oth_list->{$a1}->{$s[$c]} > $cut_off_sp_val{$a1})
					{ 
						if ($show_cnt)
						{
							$tmp .= $s[$c] . "(" . $oth_list->{$a1}->{$s[$c]} . "," .$oth_list_cnt->{$a1}->{$s[$c]}.");"; $c++;
						}
						else
						{
							$tmp .= $s[$c] . $oth_list->{$a1}->{$s[$c]} .";"; $c++;
				
						}
					}
				}
				#print STDERR $s[0] , "(" , $oth_list->{$a1}->{$s[0]} , "\n";
			
				my $nameID = $out_id;
				if ($new_names{$nameID}) { $nameID = $new_names{$nameID}; }
	
				my $qnew = $nameID; if ($q_tax) { if ($q_tax->{$a1}) { $qnew .= "\t". $q_tax->{$a1}; } else {$qnew .= "\tNA";}}
				my $bad_sec;
				if ($best2_list{$a1}->[0]->{id} - $best2_list{$a1}->[1]->{id} < $second_diff) 
				{
					if ($verb)  
					{
						print STDERR "$a1\t", $best2_list{$a1}->[0]->{id}, "\t", $best2_list{$a1}->[1]->{id}, "\t", $name{$best2_list{$a1}->[0]->{ref}}->{name}, "\t", $name{$best2_list{$a1}->[1]->{ref}}->{name}, "\n";
					}
					$bad_sec = 1;
				}
				print $fo "$qnew\t", $best_id->{$a1}->{name}, "\t";
				if ($prev) { print $fo $prev_val{$a1}, "\t"; }
				if (!$bad_sec)
				{
					if ($show_cnt)
					{
						print $fo $name{$best_id->{$a1}->{name}}->{name}, "\t", clean_name($name{$best_id->{$a1}->{name}}->{name}, $best_id->{$a1}->{id}, $cut_off), "\t(", $best_id->{$a1}->{id}, ",",$best_id->{$a1}->{cnt}, ")\t", $cut_off_sp_val{$a1},"\t", ($c-1), "\t$tmp\t", $s[$c],"\t", $oth_list->{$a1}->{$s[$c]} ,"\n";
					}
					else
					{
						print $fo $name{$best_id->{$a1}->{name}}->{name}, "\t", clean_name($name{$best_id->{$a1}->{name}}->{name}, $best_id->{$a1}->{id}, $cut_off), "\t", $best_id->{$a1}->{id}, "\t", $cut_off_sp_val{$a1},"\t", ($c-1), "\t$tmp\t", $s[$c],"\t", $oth_list->{$a1}->{$s[$c]} ,"\n";
					}
				}
				else
				{
					if ($show_cnt)
					{
						print $fo "No Clean Hit\tNoClean Hit\t", $cut_off_sp_val{$a1},"\t", ($c-1), "\t", clean_name($name{$best_id->{$a1}->{name}}->{name}, $best_id->{$a1}->{id}, $cut_off), "(", $best_id->{$a1}->{id}, ",",$best_id->{$a1}->{cnt}, ");$tmp\t", $s[$c],"\t", $oth_list->{$a1}->{$s[$c]} ,"\n";
					}
					else
					{
						print $fo "No Clean Hit\tNoClean Hit\t",$cut_off_sp_val{$a1},"\t", ($c-1), "\t", clean_name($name{$best_id->{$a1}->{name}}->{name}, $best_id->{$a1}->{id}, $cut_off), $best_id->{$a1}->{id}, ";$tmp\t", $s[$c],"\t", $oth_list->{$a1}->{$s[$c]} ,"\n";
					}
				
				}
			}
	
		}
	}
}

if ($temp_mash_id)
{
	`rm $temp_mash_id.msh`;
}

if ($run_phy)
{
	my $out1 = $out . "mat";
	my $out2 = $out . "tree";
	print STDERR "making tree on $out1 to $out2\n";

	`/usr/local/devel/ANNOTATION/tclarke/bin/make_phylogeny.R -i $out1 -m $run_phy -o $out2 -d 100`;
}

sub clean_name
{
    my $name = $_[0];
    my $id = $_[1];
	$name =~ s/([\[\]])//g;
	my @cut_off = split ",", $_[2];
    my @ids = (); while ($name =~ /(\S+)/g) { push @ids, $1; }
    my $ret = $ids[0];
    if ($id < 85) { return("Unclassified"); }
    if (scalar(@ids) > 1 && $id >= $cut_off[1]) { $ret .= " " . $ids[1]; }
		    else { $ret .= " sp."; }
    if (scalar(@ids) >2)
    {
	if ($ids[2] eq "complex") { $ret = $name;}
	if ($ids[2] =~ /\A(ssp|subsp|str.)([.]*)/ && $id >= $cut_off[0])
	{
	    return($ret . " ". $ids[2] . " " . $ids[3]);
	}
    }
    if ($ids[1] eq "sp" || $ids[1] eq "sp.")
    {
	$ret = $ids[0] . " sp.";
    }
	$ret =~ s/_/ /g;
    return $ret;
}

sub read_log_file($$)
{
	open(my $fh, "<", $_[0]) or die($_[0]);
	my $mat;
	while(my $v = <$fh>)
	{
		chomp($v);
		if ($v =~ /alnLen(\s)in(\s)(\S+):(\s+)([0-9\.]+),(\s)alnLen(\s)in(\s)(\S+):(\s+)([0-9\.]+)/)
		{
			my $id1 = $3; my $alnLen1 = $5; my $id2 =$9; my $alnLen2 = $11;
			$v = <$fh>; chomp($v);
			$v =~ /alnNuc(\s)in(\s)(\S+):(\s+)([0-9\.]+),(\s)alnNuc(\s)in(\s)(\S+):(\s+)([0-9\.]+)/;
	
			my $id1 = $3; 
			my $alnNuc1 = $5; 
			my $id2 = $9; 
			my $alnNuc2 = $11;
			if ($alnLen1) 
			{ 
				$mat->{1} = $alnNuc1/$alnLen1; 
			} 
			else 
			{
				$mat->{1}= 0;   
			}
			if ($alnLen2) 
			{ 
				$mat->{2} = $alnNuc2/$alnLen2; 
			} 
			else 
			{  
				$mat->{2} = 0;
			}
	
		}
	}
	return(sprintf("%.2f", 100*$mat->{$_[1]}));
}
