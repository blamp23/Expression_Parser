#!/usr/bin/env perl
# Script:  expressionParser.pl
# -------------------------------------------------------------------------------------
# This is the initial code that converts rFBA rules to R-matrix equations
# Example:
#	- Input: b0001 (NOT(ArcA OR Fnr))
#	- Output: -1 NOT_ArcA -1 NOT_Fnr +1 b0001
#
# Approach:
#	1. Generate prefix rFBA rule from provided infix (w/ paren's) notation
#	2. Generate an expression tree that structures the logical rule from the prefix 
#	   notation
#	3. Evaluate the expression tree yielding a fully evaluated rFBA rule (i.e. No
#	   paren's, no NOT's, only metabs, TFs, NOT_metabs, NOT_TFs, AND's, OR's
#	4. Parse fully evaluated rules into R-matrix equations (see above for example)
#
# Contains subroutines for tree traversal, expression tree evaluation, rFBA rule
# reformatting, rule parentheses checking, rFBA to R-matrix equation conversion and 
# infix to prefix notation conversion
#
# Algorithms (infix to prefix notation - in infix2prefix - as well as prefix to 
# expression tree - in main while loop as well as tree evaluator - in evaltree) 
# derived from:http://www.codepedia.com/1/Art_Expressions_p1
#
# by E.P. Gianchandani, A.R. Joyce, B.O. Palsson, and J.A. Papin
# Last revised:  February 25, 2009
# 
# Included as Protocol S1 for:
# Gianchandani, E.P., A.R. Joyce, B.O. Palsson, and J.A. Papin.  Functional 
# states of the genome-scale Escherichia coli transcriptional regulatory
# system.
# -------------------------------------------------------------------------------------

use FileHandle;
use Getopt::Long;

# CLI Input Handling
# Supports:
#   --rule "Gene (BOOLEAN EXPRESSION)"
#   --file path/to/file_containing_rule.txt
#   Or pass the rule directly as unflagged arguments
#     ./expressionParser_bl.pl Gene '(A AND B)'


my ($rule, $rule_file, $help);
GetOptions(
    'rule=s' => \$rule,
    'file=s' => \$rule_file,
    'help|h' => \$help,
);

if ($help) {
    print "Usage:\n";
    print "  $0 --rule \"Gene (A AND (B OR NOT(C)))\"\n";
    print "  $0 --file rule.txt\n";
    print "  $0 Gene '(A AND (B OR NOT(C)))'\n";
    exit 0;
}

if (!defined $rule) {
    if (defined $rule_file) {
        open my $fh, '<', $rule_file or die;
        local $/ = undef; # slurp
        $rule = <$fh>;
        close $fh;
    } elsif (@ARGV) {
        # Join remaining args into a single rule string to avoid strict quoting needs
        $rule = join(' ', @ARGV);
    } else {
        die 
    }
}



#----------------------------

#Parse apart gene/regulator from rule
#------------------------------------
my $gene = "";
my $fixedgene = "";
if ($rule =~ m/^\"/) {
	$rule =~ m/^(\".*?\")\s+/;
	$gene = $1;
	$fixedgene = $gene;
	$fixedgene =~ s/\s+/_/g;
	$fixedgene =~ s/\"//g;
	$rule =~ s/^$gene/$fixedgene/g;
} elsif ($rule =~ m/^(b\d{4,4})/) {
	#$rule =~ m/b{4,4}/;
	$gene = $1;
	$fixedgene = $1;
} elsif ($rule =~ /^[A-Z]/) {
	$rule =~ /^([A-Z].*?)\s+/;
	$gene = $1;
	#$fixedgene = $gene;
	#$fixedgene =~ s/\s+/_/g;
	#$fixedgene =~ s/\"//g;
	#$rule =~ s/^$gene/$fixedgene/g;
	($fixedgene = $gene) =~ s/\s+/_/g;
	$fixedgene =~ s/\"//g;
	$rule =~ s/$gene/$fixedgene/g;
} else {
	die "Problem with $gene: Format of regulator or target gene not recognized.\n";
}
print "PULLED GENE: $fixedgene\n";
#my @temp = split " ", $rule;
$rule =~ s/^$fixedgene//;
print "ASSOCIATED RULE: $rule\n";
#---------------------------------------

#Check rule for special words
#---------------------------------
&specialwordcheck($rule) or die "Problem with $rule: Failed special word check.";
#---------------------------------

#Check rule has matching paren's
#-------------------------------
&parencheck($rule) or die "Problem with $rule: Failed paren check.";
print "REFORMATTING: ".$rule."\n";
$rule = &reformatrule($rule);
print "DONE: ".$rule."\n";
#die;
&parencheck($rule) or die "Problem with $rule: Failed paren check.";
#--------------------------------

#Convert Infix to prefix rule notation
#-------------------------------------
my $prefix = &infix2prefix($rule);
print "PREFIX:".$prefix."\n";
#-------------------------------------

#Generate expression tree
my @prefix_rule = split /\s+/, $prefix;
@prefix_rule = reverse(@prefix_rule); #Step 1 from 
								#http://www.codepedia.com/1/Art_Expressions_p1 Algo
my %special = ();
$special{"NOT"} = 1;
$special{"AND"} = 1;
$special{"OR"} = 1;

my @treestack = ();
while (scalar(@prefix_rule) > 0) {
	my $token = shift @prefix_rule; #Step 2 from 
								#http://www.codepedia.com/1/Art_Expressions_p1 Algo
	#print "TOKEN: ".$token."\n";
	my %node;
	$node{DATA} = $token;  #Step 3 & 4.1 from 
						   #http://www.codepedia.com/1/Art_Expressions_p1 Algo
	$node{LEFT} = undef;
	$node{RIGHT} = undef;
	if (exists($special{$token})) { #Step 4 from 
							#http://www.codepedia.com/1/Art_Expressions_p1 Algo
		if ($token eq "AND" or $token eq "OR") {
			$node{LEFT} = pop @treestack;
			$node{RIGHT} = pop @treestack;
		} elsif ($token eq "NOT") {
			$node{LEFT} = pop @treestack;
		}
	}
	#$tree .= %node;
	push @treestack, \%node;
	#print "STACK: ";
	#foreach my $addy (@treestack) {
	#	print ${$addy}{DATA}." ";
	#} 
	#print "\n";
}
my $root = pop @treestack; #Step 6 from 
							#from http://www.codepedia.com/1/Art_Expressions_p1 Algo

#--------------------------------------

#Tree Traversals
#---------------
print "IN ORDER: ";
&in_order($root);
print "\n";

print "PRE ORDER: ";
&pre_order($root);
print "\n";

print "POST ORDER: ";
&post_order($root);
print "\n";
#------------------

#Evaluate the tree
print "EVAL TREE:\n";
print "INPUT RULE: ".$rule."\n";
my $newrule = &evaltree($root);
print "OUTPUT RULE: ".$newrule."\n";
#------------------

#Generate R-matrix equations
#----------------------------
my $eqnsref = &rule2Requation($fixedgene, $newrule);
my @eqns = @{$eqnsref};
my $rulecounter = 1;
print "PRINTING RULES:\n";
foreach my $r (@eqns) {
	print $fixedgene."_".$rulecounter.": ";
	print $r."\n";
	$rulecounter++;
}
#-----------------------------
#_________END MAIN______________________#

#Convert logical rules (fully evaluated) into R matrix-style equations
sub rule2Requation {
	my $bnum = $_[0]; #bnum/regulator corresponding to the rule
	my $rule = $_[1]; #fully evaluated rule

	my @eqs = split "OR", $rule;
	my @outeqs = ();
	foreach my $eqn (@eqs) {
		my @ops = split " ", $eqn;
		my $outeqn = "";
		for (my $i = 0; $i < scalar(@ops) ;$i++) {
			if ($ops[$i] ne "AND") {
				$outeqn = $outeqn." -1 ".$ops[$i];
			}
		}
		$outeqn =~ s/^\s+|\s+$//g;
		$outeqn = $outeqn." +1 ".$bnum;
		push @outeqs, $outeqn;
	}
	return \@outeqs; #return array of R-matrix-style equations
}

#Checks for special words in metabolite names other than as boolean operators 
sub specialwordcheck {
	my $rule = $_[0];
	if ($rule =~ m/\w+OR\w*|\w*OR\w+/) {
		return 0;
	} elsif ($rule =~ m/\w+AND\w*|\w*AND\w+/) {
		return 0;
	} elsif ($rule =~ m/\w+NOT\w*|\w*NOT\w+/) {
		return 0;
	} else {
		return 1;
	}
}

#Checks that open paren's match the same number of close paren's 
sub parencheck {
	my $rule = $_[0];
	#simple hack, use reg expressions to pull open paren's and close paren's
	my @openparens = $rule =~ /\(/g;
	#print "Open: ".scalar(@openparens)."\n";
	my @closeparens = $rule =~ /\)/g;
	#print "Close: ".scalar(@closeparens)."\n";
	return (scalar(@openparens) eq scalar(@closeparens)); 
}

#Evaluate a boolean expression tree (generated in MAIN)
sub evaltree {
	my $tree = $_[0]; #root of tree
	return unless($tree);
	if (${$tree}{DATA} eq "NOT") { #Handle NOT evaluations
		my $operand = &evaltree(${$tree}{LEFT});
		my @ops = split " ", $operand;
		my $new_operand = "";
		for ($i = 0;$i < scalar(@ops) ;$i++) {
			if ($ops[$i] eq "AND") { #Flip AND's to OR's (De Morgan's Law)
				$new_operand = $new_operand." "."OR";
			} elsif ($ops[$i] eq "OR") { #Flip OR's to AND's (De Morgan's Law)
				$new_operand = $new_operand." "."AND";
			} elsif ($ops[$i] =~ /^NOT_/) { #Reverse NOT_<metabs> or NOT_<regs>
				my $o = $ops[$i];
				$o =~ s/^NOT_//g;
				$new_operand = $new_operand." ".$o;
			} else {
				$new_operand = $new_operand." "."NOT_".$ops[$i]; #Generate 
													#NOT_<metabs> or NOT_<regs>
			}
		}
		$new_operand =~ s/^\s+|\s+$//g;
		print $new_operand."\n";
		return $new_operand;

	} elsif (${$tree}{DATA} eq "AND") { #Evaluate AND cases
		my $left_operand = &evaltree(${$tree}{LEFT});
		my $right_operand = &evaltree(${$tree}{RIGHT});
		print "AND: left = $left_operand | right = $right_operand\n";
		my $new_operand = "";
		#OR's are a special case that need to be fully evaluated in order to 
		#retain original intent (similar to distributive property in algebra)
		if ($left_operand =~ /\sOR\s/ or $right_operand =~ /\sOR\s/) { 
			my @lors = split " OR ", $left_operand;					   
			my @rors = split " OR ", $right_operand;				   
			my %newrules = ();										   
			for (my $i=0;$i < scalar(@lors) ;$i++) {				   
				for (my $j=0; $j < scalar(@rors) ;$j++) {
					my $newr = $lors[$i]." AND ".$rors[$j];
					$newrules{$newr} = 1;;
				}
			}
			$new_operand = join(" OR ", keys %newrules);
		} else {
			$new_operand = $left_operand." AND ".$right_operand;
		}
		print $new_operand."\n";
		return $new_operand;
	} elsif (${$tree}{DATA} eq "OR") { #OR's are easy; just glob left subtree 
									   #to right subtree seperated by an OR
		my $left_operand = &evaltree(${$tree}{LEFT});
		my $right_operand = &evaltree(${$tree}{RIGHT});
		my $new_operand = $left_operand." OR ".$right_operand;
		#print $new_operand."\n";
		return $new_operand;
	} else { #just return the value of operands
		print ${$tree}{DATA}."\n";
		return ${$tree}{DATA};
	}
}

#IN ORDER tree traversal
#Visit Left, Print Value, Visit Right
sub in_order {
	my $tree = $_[0]; #root of tree
	return unless $tree;
	&in_order(${$tree}{LEFT});
	print ${$tree}{DATA}, " ";
	&in_order(${$tree}{RIGHT});
}

#PRE ORDER tree traversal
#Print Value, Visit Left, Visit Right
sub pre_order {
	my $tree = $_[0]; #root of tree
	return unless $tree;
	print ${$tree}{DATA}, " ";
	&pre_order(${$tree}{LEFT});
	&pre_order(${$tree}{RIGHT});
}

#POST ORDER tree traversal
#Visit Left, Visit Right, Print Value 
sub post_order {
	my $tree = $_[0]; #root of tree
	return unless $tree;
	&post_order(${$tree}{LEFT});
	&post_order(${$tree}{RIGHT});
	print ${$tree}{DATA}, " ";
}

#Convert infix rules (w/ paren's) to prefix rule definition
sub infix2prefix {
	my $rule = $_[0]; #reformatted rFBA style rule (infix with parentheses, 
					  #all components delimited by a space)
	
	my @splitrule = split /\s+/, $rule;	
	my @revsplitrule = reverse(@splitrule); #step 1 from 
						#http://www.codepedia.com/1/Art_Expressions_p1 Algo
	
	my %special = ();
	$special{"NOT"} = 1;
	$special{"("} = 1;
	$special{")"} = 1;
	$special{"AND"} = 1;
	$special{"OR"} = 1;

	my $prefix = ""; #output string
	my @opstack = (); #operator stack
	while (scalar(@revsplitrule) > 0) {
		my $item = shift @revsplitrule; #step 2 from 
						#http://www.codepedia.com/1/Art_Expressions_p1 Algo 
		#print $item."\n";
		if (!exists($special{$item})) { #step 3 from 
						#http://www.codepedia.com/1/Art_Expressions_p1 Algo
			$prefix = $prefix." ".$item;
		} else {
			if ($item eq ")") { #step 4 from 
						#http://www.codepedia.com/1/Art_Expressions_p1 Algo
				push @opstack, $item;
			}
			if ($item eq "(") { #steo 6 from 
						#http://www.codepedia.com/1/Art_Expressions_p1 Algo
				my $op = pop @opstack;
				while ($op ne ")") {
					$prefix = $prefix." ".$op;
					$op = pop @opstack;
				}
			}
			# Do rules of precedence here! infinite loop implementing step 5 of algo
			# from http://www.codepedia.com/1/Art_Expressions_p1
			if ($item eq "NOT" or $item eq "AND" or $item eq "OR") {
				for (; ;) {
					if ($item eq "NOT") {
						push @opstack, $item;
						last;
					} elsif ($item eq "AND") {
						if ($opstack[scalar(@opstack)-1] eq "NOT") {
							$op = pop @opstack;
							$prefix = $prefix." ".$op;
						} else {
							push @opstack, $item;
							last;
						}
					} elsif ($item eq "OR") {
						if ($opstack[scalar(@opstack)-1] eq "NOT" 
								or $opstack[scalar(@opstack)-1] eq "AND") {
							$op = pop @opstack;
							$prefix = $prefix." ".$op;
						} else {
							push @opstack, $item;
							last;
						}
					}
				}
			}
		} 
	}
	my $op = pop @opstack;
	while ($op ne ")" and scalar(@opstack) > 0) { #step 8 from 
								#http://www.codepedia.com/1/Art_Expressions_p1 Algo
		$prefix = $prefix." ".$op;
		$op = pop @opstack;
	}
	my @prefixarray = split /\s+/, $prefix;
	$prefix = join(" ",reverse(@prefixarray)); #step 9 from 
								#http://www.codepedia.com/1/Art_Expressions_p1 Algo
	return $prefix;
}

#reformat rFBA rules such that every component is seperated by a space in order 
#to facilitate the infix2prefix subroutine
sub reformatrule {
	my $rule = $_[0]; #normal rFBA rule (infix w/ paren's)

	my @special = ();
	$special[0] = '\(';
	$special[1] = '\)';
	$special[2] = 'NOT';
	$special[3] = 'AND';
	$special[4] = 'OR';
	
	my @quotedregs = $rule =~ m/(\".*?\")/g;
	foreach my $qreg (@quotedregs) {
		print "Quoted Reg: $qreg\n";
		my $fixedqreg = $qreg;
		$fixedqreg =~ s/\s+/_/g;
		$fixedqreg =~ s/\"//g;
		print "Fixed Quoted Reg: $fixedqreg\n";
		$rule =~ s/$qreg/$fixedqreg/g;
	}
	$rule =~ s/\s+//g; #pre-process rule by removing all space characters first
	$rule =~ s/>0/_gt_0/g; #also redefine "> 0" variables 
	$rule =~ s/<0/_lt_0/g; #and "< 0" variables
	
	#Now loop through special characters introducing spaces between all components 
	#on each split
	foreach my $s (@special) {
		my @list = split /$s/, $rule,-1;
		#print "$s (".scalar(@list)."): ".join(" | ", @list)."\n";
		my $newrule = $list[0];
		for (my $i = 1; $i < scalar(@list) ;$i++) {
			if ($s eq '\(') {
				$newrule = $newrule." ( ".$list[$i];
			} elsif ($s eq '\)') {
				$newrule = $newrule." ) ".$list[$i];
			} elsif ($s  eq "NOT") {
				if ($list[$i] =~ /^_/) { #handle tricky NOT_<metab> or NOT_<TF> case
					$newrule = $newrule." ".$s.$list[$i];
				} else {
					$newrule = $newrule." ".$s." ".$list[$i];
				}
			} else {
				$newrule = $newrule." ".$s." ".$list[$i];
			}
		}
		$rule = $newrule;
		#print $rule."\n";
	}
	return $rule; #return reformatted rule
}

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 ARGUMENTS

=head1 EXAMPLE

=head1 NOTES

=head1 AUTHOR

Andrew R. Joyce (ajoyce@ucsd.edu)

=head1 COPYRIGHT

This program is free for any purpose, diabolical or otherwise.
Developed in September 2006.

=cut
