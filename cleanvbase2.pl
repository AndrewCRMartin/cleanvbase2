#!/usr/bin/perl -s
use strict;
use fasta;

########################################################################
%::usedIDs     = ();
$::FINDFIRST   = 1;
$::NOFINDFIRST = 0;

########################################################################
UsageDie() if(defined($::h));

if(open(my $fp, '<', $ARGV[0]))
{
    my %infos = ();
    my %seqs  = ();

    while(my($id, $info, $seq) = fasta::ReadFasta($fp))
    {
        last if ($id eq '');
        $id =~ s/^\s+//;
        $infos{$id} = $info;
        $seqs{$id}  = $seq;
    }
    close($fp);

    # Extract the canonical names from the IDs
    my %keyIDs = BuildKeyIDs(keys %infos);

    if(defined($::genes) || defined($::noalleles))
    {
        %keyIDs = RemoveNonGenes(%keyIDs);
    }

    if(defined($::noalleles))
    {
        %keyIDs = RemoveAlleles(%keyIDs);
    }


    foreach my $id (sort keys %keyIDs)
    {
        my $printID = $id;
        $printID =~ s/\>\s+//;
        my $printInfo = $infos{$id};
        $printInfo =~ s/\>\s+//;
        if($printInfo eq $printID)
        {
            $printInfo = '';
        }
        else
        {
            $printInfo = "|$printInfo";
        }
        print ">$keyIDs{$id}|$printID$printInfo\n";
        print "$seqs{$id}\n" if(!defined($::noseqs));
    }

}


########################################################################
# %keyIDs = BuildKeyIDs(@ids)
# ---------------------------
# Takes an array of IDs such as
# '> humIGHV057 IGHV3-66*01, IGHV3-66*04, 8-1B+, MTGLa 291 bp'
# and builds a hash of the best short ID (in this case 'IGHV3-66*01')
# indexed by the full ID.
#
# The routine performs several passes first trying to find proper
# allele gene names and then things containing IG[HKL]V followed by
# DP and \d-\d. The first field in the ID is not searched but used
# as a last resort if none of the matches worked.
#
# 11.10.16 Original   By: ACRM
sub BuildKeyIDs
{
    my(@ids) = @_;

    my %keyIDs  = ();

    foreach my $id (@ids)
    {
        my @fields   = split(/[\,\s+]+/, $id);

        # Pass 1, look for a field with a *
        if(!LookFor('\\*', $id, \@fields, \%keyIDs, $::FINDFIRST))
        {
            # Pass 2, look for a field containing IG.V
            if(!LookFor('IG[HKL]V', $id, \@fields, \%keyIDs, $::NOFINDFIRST))
            {
                # Pass 3, look for a field containing DP
                if(!LookFor('DP', $id, \@fields, \%keyIDs, $::NOFINDFIRST))
                {
                    # Pass 4, look for a field containing n-n
                    if(!LookFor('\\d\-\\d', $id, \@fields, \%keyIDs, $::NOFINDFIRST))
                    {
                        # Still nothing, so just use first field
                        $keyIDs{$id} = $fields[0];
                    }
                }
            }
        }
    }

    return(%keyIDs);
}


########################################################################
# $found = LookFor('\\*', $id, \@fields, \%keyIDs, $::FINDFIRST)
# --------------------------------------------------------------
# Takes the $id and the $id pre-split into fields ($aFields)
# Searches each of the fields (skipping the first) to see if it matches 
# $regex. 
#
# If $findFirst is not set, then it simply stores the field, indexed
# by the $id, in the $hKeyIDs hash. It also records that this field
# has been used by storing it in the global %::usedIDs hash so that if
# this is encountered again, the first field will be appended onto
# this field. For example if DP-61 is found as a matched field twice,
# then the first one will appear simply as DP-61, the second as
# DP-61*humIGHV085
# 
# If $findFirst is set, then rather than simply taking the first field
# that matches, all fields that match are checked and the one that is
# alphabetically first is selected. For example the $id
# '> humIGHV057 IGHV3-66*01, IGHV3-66*04, 8-1B+, MTGLa 291 bp' when
# using the $regex /*/ would match IGHV3-66*01 and IGHV3-66*04.
# This would ensure that IGHV3-66*01 was selected rather than
# IGHV3-66*04.
#
# This routine is used in several passes as a hierarchy for assigning
# key IDs for the entries. Fields with a * are proper genes with allele
# numbers so would be searched for first.
#
# 11.10.16 Original   By: ACRM
sub LookFor
{
    my($regex, $id, $aFields, $hKeyIDs, $findFirst) = @_;
    my $field    = '';
    my $fieldNum = 0;

    if($findFirst)
    {
        my @matches = ();

        foreach my $field (@$aFields)
        {
            if($fieldNum)
            {
                if($field =~ /$regex/)
                {
                    $field =~ s/\.\.\.$//;
                    push @matches, $field;
                }
            }
            $fieldNum++;
        }

        if(scalar(@matches))
        {
            my $firstMatch = "ZZZZZZZZZZZZZZ";
            foreach my $match (@matches)
            {
                if($match lt $firstMatch)
                {
                    $firstMatch = $match;
                }
            }

            if(defined($::usedIDs{$firstMatch}))
            {
                $firstMatch .= "*" . $$aFields[0];
            }
        
            $$hKeyIDs{$id}          = $firstMatch;
            $::usedIDs{$firstMatch} = $id;
            return(1);
        }
    }
    else
    {
        foreach my $field (@$aFields)
        {
            if($fieldNum)
            {
                if($field =~ /$regex/)
                {
                    $field =~ s/\.\.\.$//;
                    if(defined($::usedIDs{$field}))
                    {
                        $field .= "*" . $$aFields[0];
                    }
        
                    $$hKeyIDs{$id}     = $field;
                    $::usedIDs{$field} = $id;
                    return(1);
                }
            }
            $fieldNum++;
        }
    }
    return(0);
}



########################################################################
# %keyIDs = RemoveNonGenes(%keyIDs)
# ---------------------------------
# Takes a hash of key IDs (e.g. IGHV4-31*02, humIGHV309, DP-24 or 
# V3-22P) indexed by the full ID
# (e.g. '> humIGHV197 IGHV4-31*02, VH4-31, V4-31+ 297 bp') and removes
# all the entries that do not have proper gene names like IGHV4-31*02
#
# This is designed to remove all pseudogenes, etc. from the returned
# hash
#
# 11.10.16 Original   By: ACRM
sub RemoveNonGenes
{
    my(%keyIDs) = @_;
    my %keptIDs = ();

    foreach my $id (keys %keyIDs)
    {
        if($keyIDs{$id} =~ /^IG[HKL]V/)
        {
            $keptIDs{$id} = $keyIDs{$id};
        }
    }

    return(%keptIDs);
}


########################################################################
# %keyIDs = RemoveAlleles(%keyIDs)
# --------------------------------
# Takes a hash of key IDs (e.g. IGHV4-31*02) indexed by the full ID
# (e.g. '> humIGHV197 IGHV4-31*02, VH4-31, V4-31+ 297 bp') and removes
# allelic variants (e.g. if both IGHV4-31*02 and IGHV4-31*03 were in
# the hash, the returned hash would only contain IGHV4-31*02. The
# allelic vaiant with the lowest suffix number is retained.
#
# 11.10.16 Original   By: ACRM
sub RemoveAlleles
{
    my(%keyIDs)  = @_;
    my %bestKeyIDs  = ();
    my %suffixes = ();

    foreach my $id (keys %keyIDs)
    {
        my @fields = split(/\*/, $keyIDs{$id});
        my $stem   = $fields[0];
        my $suffix = $fields[1];

        if(defined($suffixes{$stem}))
        {
            if($suffix lt $suffixes{$stem})
            {
                $suffixes{$stem} = $suffix;
                DeleteEntry(\%bestKeyIDs, $stem);
                $bestKeyIDs{$id} = $keyIDs{$id};
            }
        }
        else
        {
            $bestKeyIDs{$id} = $keyIDs{$id};
            $suffixes{$stem} = $suffix;
        }
    }

    return(%bestKeyIDs);
}

########################################################################
# DeleteEntry(\%bestKeyID, $stem)
# -------------------------------
# Takes a reference to a hash indexed by ID and containing the allele
# name (e.g. IGHV4-31*02) together with a stem (e.g. IGHV4-31).
# The code searches the hash for a value starting with the spcified
# stem and deletes this entry from the hash.
# This is used so that the entry can be replaced by a different
# allele number (e.g. IGHV4-41*01)
#
# 11.10.16 Original   By: ACRM
sub DeleteEntry
{
    my($hBestKeyIDs, $keyID) = @_;
    foreach my $id (keys %$hBestKeyIDs)
    {
        if($$hBestKeyIDs{$id} =~ /^$keyID\*/)
        {
            delete $$hBestKeyIDs{$id};
            last;
        }
    }
}

########################################################################
sub UsageDie
{
    print <<__EOF;

cleanvbase2 V1.0 (c) 2016, UCL, Dr. Andrew C.R. Martin

Usage: cleanvbase2 [-genes][-noalleles][-noseqs] vbase2.faa
       -genes     Only display proper genes identified as IGxVn-m*aa
       -noalleles Keep only the lowest numbered allelic variant of
                  proper genes
       -noseqs    Show only the FASTA headers, not the actual sequences

Cleans up a VBASE2 Fasta file to produce a file with identifiers using
standard gene names instead of VBASE2's own identifiers.

Also allows only proper genes to be selected (discarding pseudogenes, 
etc.) and also only one example of each gene where there are allelic
variants.

__EOF

    exit 0;
}
