#!/usr/bin/awk -f


#define NEUTRAL 0
#define NEUTRAL_ACCEPTOR 1
#define NEUTRAL_DONOR 2
#define NEUTRAL_ACCEPTOR_DONOR 3
#define POSITIVE_ACCEPTOR 4
#define POSITIVE_DONOR 5
#define POSITIVE_ACCEPTOR_DONOR 6
#define NEGATIVE_ACCEPTOR 7
#define NEGATIVE_DONOR 8
#define NEGATIVE_ACCEPTOR_DONOR 9

$4==0&&$5==0&&$6==0{print 0}
$4==1&&$5==0&&$6==0{print 1}
$4==0&&$5==1&&$6==0{print 2}
$4==1&&$5==1&&$6==0{print 3}

$4==1&&$5==0&&$6==1{print 4}
$4==0&&$5==1&&$6==1{print 5}
$4==1&&$5==1&&$6==1{print 6}

$4==1&&$5==0&&$6==-1{print 7}
$4==0&&$5==1&&$6==-1{print 8}
$4==1&&$5==1&&$6==-1{print 9}
