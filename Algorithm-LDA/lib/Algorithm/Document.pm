#!/usr/bin/perl
#
# @File Document.pm
# @Author prometheus
# @Created Jul 5, 2017 10:54:51 AM
#

package Algorithm::Document;

use strict;
use 5.006;
use warnings FATAL => 'all';

use base qw(Class::Accessor::Fast);
Algorithm::Document->mk_accessors("_length");
Algorithm::Document->mk_accessors("_total");
Algorithm::Document->mk_accessors("_id");
Algorithm::Document->mk_accessors("_words");
Algorithm::Document->mk_accessors("_counts");
Algorithm::Document->mk_accessors("_topics");

my $self;

sub new 
{
    my $class = shift;
    $self = {
       
       _length => 0,
       _total => 0,
       _id => 0,
       _words => [],
       _counts =>[],
       _topics => {},
    };

    
    bless $self, $class;
            
    return $self;
}

1;