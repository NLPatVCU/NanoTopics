#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Perl::LDA' ) || print "Bail out!\n";
}

diag( "Testing Perl::LDA $Perl::LDA::VERSION, Perl $], $^X" );
