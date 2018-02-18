#!/usr/bin/perl
#
# @File LDA.pm
# @Author prometheus
# @Created Jul 12, 2017 3:28:49 PM
#

package Algorithm::LDA;

use strict;
use 5.006;
use warnings FATAL => 'all';
use Algorithm::Document;

use List::Util qw(sum);

use constant pi => 4*atan2(1, 1);
use constant e  => exp(1);



my $initialAlpha;
my $K;
my $k;
my $maxiter;
my $initMethod;
my $convFactor;

my $alpha;

my $data;
my $vocab;
my $output;

my @docs;
my $nterms;
my $ndocs;

my %topicMap = ();
my %documentTopics = ();
my %documentMap = ();
my %topicWords = ();

my %topics = ();

my %vocabulary = ();

my @phi = ();
my @theta = ();
my @beta = ();

my $self;

sub new 
{
    my $class = shift;
    $self = {
        _alpha => shift,
        _ntopics => shift,
        _maxiter => shift,
        _initMethod => shift,
        _convergence => shift,
        _data => shift,
        _vocab => shift,
        _output => shift,
    };
    
    
    $initialAlpha = $self->{_alpha};
    $K = $self->{_ntopics};
    $k = $K - 1;
    $maxiter = $self->{_maxiter};
    $initMethod = $self->{_initMethod};
    $convFactor = $self->{_convergence};
    
    
    $data = $self->{_data};
    $vocab = $self->{_vocab};
    $output = $self->{_output};
    
    bless $self, $class;
            
    return $self;
}

sub init
{      
    readData($data);
    readVocabulary($vocab);
    
    beta();
    
    for my $doc(0..$ndocs-1)
    {
        my @array = ();
        for my $i (0..$k)
        {
            push(@array, 0);
        }
        push(@theta, \@array);
    }
    
    for my $i(0..$k)
    {
        my @array = ();
        for my $word (0..$nterms - 1)
        {
            push(@array, 0);
        }
        push(@phi, \@array);
    }
    
    mkdir "Results/";
}

sub run
{
    my $converged = 0;
    my $iter = 0;
    
    $alpha = $initialAlpha;
    
    while(($converged != 1) && (($iter <= $maxiter) || $maxiter == -1))
    {
        $iter++;
        print "...Iteration: $iter...\n";
        
        my $convcount = 0;
        my $convtotal = 0;
        
        for my $d (0..$ndocs - 1)
        {
            my $doc = $docs[$d];
            
            if (($d % 1000) == 0) {
                print"...document $d...\n";
            }
            
            for my $n (0..$doc->_length - 1)
            {
                my $word = $doc->_words->[$n];
                my $count = $doc->_counts->[$n];        
                        
                #my $prevTopic = $topics{$word};
                my $prevTopic = ${$doc->_topics}{$word};
                
                $documentMap{$doc}-=$count;
                $documentTopics{$doc}{$prevTopic}-=$count;
                $topicMap{$prevTopic}-=$count;
                $topicWords{$prevTopic}{$word}-=$count;
                
                my $dist = distribution($doc, $word);
                my $newTopic = getTopic($dist);
                
                ${$doc->_topics}{$word} = $newTopic;
                $topics{$word} = $newTopic;
                
                if($prevTopic == $newTopic)
                {
                    $convcount++;
                }
                $convtotal++;
                
                $documentMap{$doc}+=$count;
                $documentTopics{$doc}{$newTopic}+=$count;
                $topicMap{$newTopic}+=$count;
                $topicWords{$newTopic}{$word}+=$count;
            }
        }
        
        print "$convcount : $convtotal\n";
        if($convcount >= $convtotal * $convFactor)
        {
            print $convcount . "\n";
            $converged = 1;
            print "Converged on iteration $iter.\n";
        }
        
        saveResults();
    }
}

sub distribution
{
    
    my($doc, $word) = @_;
        
    my @distrib = ();
    my $dist = 0.0;
    
    for my $i (0..$k)
    {
        my $tw = (exists $topicWords{$i}{$word}) ? $topicWords{$i}{$word} : 0;
        my $t = (exists $topicMap{$i}) ? $topicMap{$i} : 0;
        my $dt = (exists $documentTopics{$doc}{$i}) ? $documentTopics{$doc}{$i} : 0;
        my $d = (exists $documentMap{$doc}) ? $documentMap{$doc} : 0;
        
        my $p = ($tw + $beta[$i][$word]) / 
                    ($t + $beta[$i][$word] * $nterms);                   
        my $th = ($dt + $alpha) / 
                        ($d + $alpha * $K);
                        
        $phi[$i][$word] = $p;
        $theta[$doc->_id][$i] = $th;
                        
        $dist += $p * $th;
        
        push(@distrib, $dist);
    }
    
    return \@distrib;
}

sub getTopic
{
    my($distrib) = @_;
    
    
    my $sample = rand($distrib->[$k]);
    
    for my $i (0..$k)
    {
        if($sample < $distrib->[$i])
        {
            return $i;
        }
    }
    
    return $k;
}

sub beta
{   
    
    my $e_value = 1.0 / $K;
            
    for my $i (0..$k) 
    {
        for my $n (0..$nterms - 1) 
        {
            $beta[$i][$n] = $e_value;
        }
    }
      
    for(1..1000000) 
    {
        my $d = rand()* $e_value;
        my $i = int(rand($K));
        my $n1 = int(rand($nterms));
        my $n2 = int(rand($nterms));
               
        $beta[$i][$n1]+=$d;
        $beta[$i][$n2]-=$d;
               
        if($beta[$i][$n2] <= 0 || $beta[$i][$n1] >=1)
        {
            $beta[$i][$n1]-=$d;
            $beta[$i][$n2]+=$d;              
        }
    }
    
#    for my $i (0..$k) 
#    {
#        my $i_star = int(rand($V))-1;
#        for my $m (0..$v) 
#        {
#            $beta[$i][$m] = $e_value / $v;
#            $beta[$i][$m]  = 1-$e_value if $m==$i_star;
#                              
#        }
#    }
}

sub saveResults
{
    my $dir = "Results/$output";
    
    my $file = "$dir.topics";
    open(my $fh, '>', $file) or die "Could not open file '$file' $!";
    
    for my $word (keys %topics)
    {
        print $fh $word . " " . $topics{$word} . " " . $phi[$topics{$word}][$word] . "\n";
    }
    close $fh;
    
    $file = "$dir.assignments";
    open(my $fh1, '>', $file) or die "Could not open file '$file' $!";
    
    #my @keys = sort {$topics{$a} <=> $topics{$b}} keys %topics;
    my %reverse = ();
    
    for my $key (keys %topics)
    {
        if(not exists $reverse{$topics{$key}})
        {
            $reverse{$topics{$key}} = ();
        }
        push(@{$reverse{$topics{$key}}}, $key);
    }
    
    for my $key (sort {$a <=> $b}keys %reverse)
    {
        print $fh1 "\n" . $key . "\n";
        my @array = sort { $phi[$key][$b] <=> $phi[$key][$a] } @{$reverse{$key}};
        for my $i (0..5)
        {
            if($i <= $#array)
            {
                print $fh1 $vocabulary{$array[$i]} . " " . $phi[$key][$array[$i]] . "\n";
            }
        }
    }
    
    close($fh1);
    
    
    my $words = "test/word-assignmentsAPSmall.dat";
    open(my $fh2, '>', $words) or die "Could not open file '$words' $!";
    
    for my $d (0..$ndocs - 1)
    {
        my $doc = $docs[$d];
        print $fh2 $doc->_length;
        for my $n (0..$doc->_length - 1)
        {
            my $word = $doc->_words->[$n];
            my $count = $doc->_counts->[$n];
            #my $topic = sprintf("%02d", $topics{$word});
            my $topic = sprintf("%02d", ${$doc->_topics}{$doc->_words->[$n]});
            print $fh2 " " . sprintf("%04d", $word) . ":" . sprintf("%04d", $count) . ":" . $topic;
        }
        print $fh2 "\n";
    }
}

sub readData
{
    my($filename) = @_;
    
    open(my $fh, $filename) or die "Could not open file '$filename' $!";
    my $length = 0;
    my $count = 0;
    my $word = 0;
    my $n = 0;
    my $nd = 0;
    my $nw = 0;
    
    @docs = ();
    $ndocs = 0;
    $nterms = 0;
    
    while (my $doc = <$fh>) {
        chomp $doc;
        my @arr = split(/\s+/, $doc);
        $length = $arr[0];
        
        $docs[$nd] = new Algorithm::Document();
        $docs[$nd]->_id($nd);
        #$docs[$nd]->set("_length", $length);
        $docs[$nd]->_length($length);

        for my $i (1..$length) {
            my @arr1 = split(/:/, $arr[$i]);
            my $word = $arr1[0];
            my $count = $arr1[1];
            push(@{$docs[$nd]->_words}, $word);
            push(@{$docs[$nd]->_counts}, $count);
            #$docs[$nd]->set("_total", $docs[$nd]->_total + $count);
            $docs[$nd]->_total($docs[$nd]->_total + $count);
	    if ($word >= $nw) {
                $nw = $word + 1;
            }
            
            my $topic = int(rand($K));
            ${$docs[$nd]->_topics}{$word} = $topic;
            #$topics{$word} = $topic;
            
            $documentMap{$docs[$nd]}+=$count;
            $documentTopics{$docs[$nd]}{$topic}+=$count;
            $topicMap{$topic}+=$count;
            $topicWords{$topic}{$word}+=$count;
        }
        $nd++;
    }
    
    $ndocs = $nd;
    $nterms = $nw;
    
    print "ndocs: $ndocs\nnterms: $nterms\n";
}

sub readVocabulary
{
    my($filename) = @_;
    
    open(my $fh, $filename) or die "Could not open file '$filename' $!";
    
    my $i = 0;
    while (my $word = <$fh>) {
        chomp $word;
        $vocabulary{$i} = $word;
        $i++;
    }
}

1;
