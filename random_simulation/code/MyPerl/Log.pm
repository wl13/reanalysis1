#!/usr/bin/perl -w
#
#   Log.pm -- Put a log about the feedback of the program
#   Author: Nowind
#   Created: 2010-11-13
#   Last Modified: 2010-11-13

use strict;

package MyPerl::Log;

use Symbol;

sub new
{
    my ($pkg, $filename) = @_;
    #my $obj = gensym();
    my $obj = \*STDERR;
    open ($obj, ">> $filename") || return undef;
    bless $obj, $pkg;
}

sub start
{
    my $obj  = shift;
    
    my $time = scalar localtime();

    print $obj "=" x 24, "[ $time ]",
               "=" x 24, "\n";
}

sub error
{
    my ($obj, $error) = @_;
    
    my $time = scalar localtime();
    
    print $obj "$time: $error";
}

sub end {
    my $obj  = shift;
    
    my $time = scalar localtime();
    
    print $obj "$time\n";
    print $obj "=" x 76 . "\n";
    
    close $obj;
}

1;
