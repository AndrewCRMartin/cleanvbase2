package fasta;
#*************************************************************************
#
#   Program:    
#   File:       fasta.pm
#   
#   Version:    V1.0
#   Date:       05.11.13
#   Function:   Read FASTA entries
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 2013
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#   while((($id, $info, $sequence) = fasta::ReadFasta(*STDIN)) && ($id ne ""))
#   {
#       print "ID: '$id'\nINFO: '$info'\nSEQ: '$sequence'\n";
#   }
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0   05.11.13   Original   By: ACRM
#
#*************************************************************************
$::RFPrevHeader = "";
#*************************************************************************
sub ReadFasta
{
    my($fp) = @_;

    my $seq = "";
    while(<$fp>)
    {
        chomp;
        s/^\s+//;

        if(length)
        {
            if(/^>/)
            {
                if($::RFPrevHeader eq "")
                {
                    $::RFPrevHeader = $_;
                }
                else
                {
                    my $id = GetFastaID($::RFPrevHeader);
                    my $info = $::RFPrevHeader;
                    $::RFPrevHeader = $_;
                    return($id, $info, $seq);
                }
            }
            else
            {
                s/\s//g;
                $seq .= $_;
            }
        }
    }

    if($::RFPrevHeader eq "")
    {
        return("","","");
    }

    my $id = GetFastaID($::RFPrevHeader);
    my $info = $::RFPrevHeader;
    $::RFPrevHeader = "";

    return($id, $info, $seq);
}

#*************************************************************************
sub PrintFasta
{
    my($fh, $header, $seq, $width) = @_;
    $width = 60 if(!defined($width) || ($width == 0));

    print $fh "$header\n";
    while(length($seq))
    {
        print $fh substr($seq, 0, $width) . "\n";
        $seq = substr($seq, $width);
    }
}

#*************************************************************************
sub GetFastaID
{
    my($text) = @_;

    $text =~ s/^>//;
    my(@parts) = split(/\|/, $text);
    if(scalar(@parts) > 1)
    {
        if(length($parts[0])<3)
        {
            return($parts[1]);
        }
        else
        {
            return($parts[0]);
        }
    }

    return($text);
}

1;
