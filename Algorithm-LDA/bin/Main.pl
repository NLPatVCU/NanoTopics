#!/usr/bin/perl
#
# @File Main.pl
# @Author prometheus
# @Created Jul 5, 2017 10:10:25 AM
#

use strict;
use warnings;
use Algorithm::LDA;


my $alpha = 0.1;
my $ntopics = 2;
my $maxiter = -1;
my $initmethod = "random";
my $convergence = 0.85;
my $data = "test/ap.dat";
my $vocab = "test/vocab.txt";
my $output = "documentsAP";

my $lda = new Algorithm::LDA($alpha, $ntopics, $maxiter, $initmethod, $convergence, $data, $vocab, $output);
$lda->init();
$lda->run();
