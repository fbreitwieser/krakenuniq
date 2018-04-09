package File::SortedSeek;
use strict;
use warnings;
use Time::Local;
require Exporter;

use vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION );

@ISA         = qw( Exporter );
@EXPORT      = ();
@EXPORT_OK   = qw( alphabetic numeric find_time get_between get_last );
%EXPORT_TAGS = ( 'all' => \@EXPORT_OK );
$VERSION     = '0.015';

my $descending  = 0;
my $cuddle      = 0;
my $line_length = 80;
my $error_msg   = '';
my $silent      = 0;
my $exact_match = 0;
my %months = ( Jan => 0, Feb => 1, Mar => 2, Apr => 3,
               May => 4, Jun => 5, Jul => 6, Aug => 7,
               Sep => 8, Oct => 9, Nov => 10,Dec => 11);
my $default_rec_sep = qw/\015\012|\015|\012/;

# some subs to set optional vars OO style
sub set_cuddle      { $cuddle = 1 };
sub set_no_cuddle   { $cuddle = 0 };
sub set_descending  { $descending = 1 };
sub set_ascending   { $descending = 0 };
sub set_max_tries   { };  # legacy method, no effect
sub set_line_length { };  # legacy method, no effect
sub set_silent      { $silent = 1 };
sub set_verbose     { $silent = 0 };
sub was_exact       { $exact_match };
sub error           { $error_msg; };

sub _alphabetic_compare { $descending ? $_[1] cmp $_[0] : $_[0] cmp $_[1] }

sub alphabetic {
    local *FILE     = shift;
    my $string      = shift;
    my $munge_ref   = shift;
    $error_msg   = '';
    _look( *FILE, $string, \&_alphabetic_compare, $munge_ref );
}

sub _numeric_compare { $descending ? $_[1] <=> $_[0] : $_[0] <=> $_[1] }

sub numeric {
    local *FILE     = shift;
    my $number      = shift;
    my $munge_ref   = shift;
    $error_msg   = '';
    _look( *FILE, $number, \&_numeric_compare, $munge_ref  );
}

sub find_time {
    local *FILE     = shift;
    my $find        = shift || time;
    my $not_gmtime  = shift;
    $error_msg   = '';
    $find = get_epoch_seconds($find,$not_gmtime) unless $find =~ m/^[\d.]+$/;
    _look( *FILE, $find, \&_numeric_compare, \&get_epoch_seconds );
}

sub get_epoch_seconds {
    my ($line, $not_gmtime) = @_;
  return undef unless defined $line;
    my ($wday,$mon,$mday,$hours,$min,$sec,$year);
    # look for asctime format: Tue May 27 15:45:00 2008
    # ignore wday token as this is often dropped ie linux kernel messages
    if ($line =~ m/(\w{3})\s+(\d{1,2})\s+(\d\d):(\d\d):(\d\d)\s+(\d{4})/ ) {
       ($mon,$mday,$hours,$min,$sec,$year) = ($1,$2,$3,$4,$5,$6);
    }
    # look for apache time format: [21/May/2008:17:49:39 +1000]
    # ignore the time offset
    elsif($line =~ m!\[(\d{1,2})/(\w{3})/(\d{4}):(\d\d):(\d\d):(\d\d)!x ) {
       ($mday,$mon,$year,$hours,$min,$sec) = ($1,$2,$3,$4,$5,$6);
    }
    # look for straight epochtime data (ie squid log)
    elsif($line =~ m/^(\d+)/) {
        return $1;
    }
    unless ($year) {
        $error_msg = "Unable to find time like string in line:\n$line";
        warn $error_msg unless $silent;
      return undef;
    }
    $mon = $months{$mon};   # convert to numerical months 0 - 11
  return $not_gmtime ? timelocal($sec,$min,$hours,$mday,$mon,$year):
                       timegm($sec,$min,$hours,$mday,$mon,$year);
}

sub get_between {
    local *FILE = shift;
    my $begin   = shift || 0;
    my $finish  = shift || 0;
    my $rec_sep = shift || $default_rec_sep;
    $error_msg   = '';
    ($begin , $finish) =  ($finish, $begin) if $begin > $finish;
    my $bytes   = $finish - $begin;
    sysseek FILE, $begin, 0;
    my $read = sysread(FILE, my $buffer, $bytes);
    if ( $read < $bytes ) {
        $error_msg  = "Short read\nWanted: $bytes Got: $read\n";
        warn $error_msg unless $silent;
      return undef;
    }
    $buffer = substr $buffer, 0, $bytes;
    my @lines = split $rec_sep, $buffer;
  return wantarray ? @lines : [ @lines ];
}

sub get_last {
    local *FILE    = shift;
    my $num_lines  = shift;
    my $rec_sep    = shift || $default_rec_sep;
    $error_msg   = '';
    my @stat = stat(FILE) or return undef;
    my($size,$blksize) = @stat[7,11];
    $blksize ||= 8192;
    # grab the first chunk back from eof at block offset
    my $pos = $size - (($size % $blksize)|| $blksize );
    my $file = '';
    my ($buf, $lines);
    for(;;) {
        $pos = 0 if $pos < 0;
        sysseek(FILE,$pos,0);
        sysread(FILE, $buf, $blksize) or last; # returns 0 at eof;
        $file = $buf.$file;
        my $lines = () = $file =~ m/$rec_sep/g;
      last if $lines > $num_lines or $pos == 0;
        $pos -= $blksize;
    }
    my @file = split /$rec_sep/, $file;
    if ( $num_lines > @file ) {
        $error_msg = "Unable to find $num_lines\n";
        warn $error_msg unless $silent;
      return wantarray ? @file : \@file;
    }
    else {
        $num_lines = $#file - $num_lines + 1;
      return wantarray ? @file[$num_lines..$#file] : [@file[$num_lines..$#file]];
   }
}

# Modified version of Perl Search::Dict's look()

sub _look {
    local *FILE = shift;
    my($key,$comp,$xfrm) = @_;
    local $_;
    return undef if not defined $key;
    my @stat = stat(FILE) or return undef;
    my($size, $blksize) = @stat[7,11];
    $blksize ||= 8192;
    # find the right block
    my($min, $max) = (0, int($size / $blksize));
    my $mid;
    while ($max - $min > 1) {
        $mid = int(($max + $min) / 2);
        seek(FILE, $mid * $blksize, 0) or return undef;
        <FILE> if $mid; # probably a partial line
        $_ = <FILE>;
        $_ = $xfrm->($_) if $xfrm;
        chomp;
        (defined($_) && $comp->($_, $key) < 0) ? $min = $mid : $max = $mid;
    }
    # find the right line
    $min *= $blksize;
    seek(FILE,$min,0) or return undef;
    <FILE> if $min; # probably a partial line
    my $prev_min = $min;
    for (;;) {
        $min = tell(FILE);
        defined($_ = <FILE>) or last;
        $_ = $xfrm->($_) if $xfrm;
        chomp;
        my $cmp = $comp->($_, $key);
        $exact_match = $cmp==0 ? 1 : 0;
        if(!$cuddle and $cmp  >= 0){
            seek(FILE,$min,0);
            return $min;
        }
        if($cuddle and $cmp > 0){
            seek(FILE,$prev_min,0);
            return $prev_min;
        }
        $prev_min = $min;
    }
    return undef;
}

1;

__END__

=pod

=for stopwords Stig refactored ta da hh mm ss dd mm yyyy recognised

=head1 NAME

File::SortedSeek - A Perl module providing fast access to large files

=head1 SYNOPSIS

  use File::SortedSeek ':all';
  open BIG, $file or die $!;

  # find a number or the first number greater in a file (ascending order)
  $tell = numeric( *BIG, $number );
  # read a line in from where we matched in the file
  $line = <BIG>;
  print "Found exact match as $line" if File::SortedSeek:was_exact();

  # find a string or the first string greater in a file (alphabetical order)
  $tell = alphabetic( *BIG, $string );
  $line = <BIG>;

  # find a date in a logfile supplying a scalar localtime type string
  $tell = find_time( *BIG, "Thu Aug 23 22:59:16 2001" );
  # or supplying GMT epoch time
  $tell = find_time( *BIG, 998571554 );
  # get all the lines after our date
  @lines = <BIG>;

  # get the lines between two logfile dates
  $begin  = find_time( *LOG, $start );
  $end    = find_time( *LOG, $finish );
  # get lines as an array
  @lines = get_between( *LOG, $begin, $end );
  # get lines as an array reference
  $lines = get_between( *LOG, $begin, $end );

  # use you own sub to munge the file line data before comparison
  $tell = numeric( *BIG, $number, \&epoch );
  $tell = alphabetic( *BIG, $string, \&munge_line );

  # use methods on files in reverse alphabetic or descending numerical order
  File::SortedSeek::set_descending();

  # for inexact matches set FH so first value read is before and second after
  File::SortedSeek::set_cuddle();

  # get last $n lines of any file as an array
  @lines = get_last( *BIG, $n )
  # or an array reference
  $lines = get_last( *BIG, $n )
  # change the input record separator from the OS default
  @lines = get_last( *BIG, $n, $rec_sep )

=head1 DESCRIPTION

File::SortedSeek provides fast access to data from large files. Three
methods numeric() alphabetic() and find_time() depend on the file data
being sorted in some way. Dictionaries are and obvious example but log files
are also sorted (by date stamp). The get_between() method can be used to get
a chunk of lines efficiently from anywhere in the file. The required position(s)
for the get_between() method are supplied by the previous methods. The
get_last() method will efficiently get the last N lines of any file, sorted
or not.

With sorted data a linear search is not required. Here is a typical linear
search

    while (<FILE>) {
        next unless /$some_cond/
        # found cond, do stuff
    }

Remember that old game where you try to guess a number between lets say 0
and say 128? Let's choose 101 and now try to guess it.

Using a linear search is the same as going 1 higher 2 higher 3 higher ...
100 higher 101 correct! Consider a geometric approach: 64 higher 96 higher
112 lower 104 lower 100 higher 102 lower - ta da must be 101! This is the
halving the difference  or binary search method and can be applied to any data
set where we can logically say higher or lower. In other words any sorted data
set can be searched like this. It is a far more efficient method - see the
SPEED section for a brief analysis.

=head2 numeric() and alphabetic() - The two basic methods

There are two basic methods - numeric() to do numeric searches and
alphabetic() that does alphabetic searches.

You call the functions like this:

    $tell = numeric( *BIG, $find );
    $tell = alphabetic( *BIG, $find );

These methods take two required arguments. *BIG is a FILEHANDLE to read from.
$find is the item you wish to find. $find must be appropriate to the function
as the numeric method will make numeric comparisons ( == < > ). Similarly the
alphabetic method makes string comparisons ( eq lt gt ). You will get strange
results if you use the wrong method just as you do if you use < when you
actually meant lt

=head3 Return values with search success and failure

The return value from the numeric() and alphabetic() methods depend on the
result of the search. If the search fails the return value is undef.
A search can succeed in two ways. If an exact match is found then the
current file position pointer is set to the beginning of the matching line.
The return value is the corresponding response from tell(). This means that
the next read from <FILEHANDLE> will return the matching line.
Subsequent reads return the following lines as expected.

Alternatively a search will succeed if a point in the file can be found such
that $find is cuddled between two adjacent lines. For example consider
searching for the number 42 in a file like this:

    ..
    36
    40  <- Before
    44  <- After
    48
    ..

The number 42 is not actually there but the search will still succeed as it
is between 40 and 44. By default the file position pointer is set to the
beginning of the line '44' so the next read from <FILEHANDLE> will return
this line. If the File::SortedSeek::set_cuddle() function is called then the
file position pointer will be set to the beginning of line '40' so that the
first two reads from <FILEHANDLE> will cuddle the in-between value in $find.

You can find out if the match was exact by checking the value returned by
File::SortedSeek::was_exact() which will be true if the match was exact.

=head3 Adding line munging to make the basic methods more useful

Both the numeric and alphabetic subs take an optional third argument.
This optional argument is a reference to a subroutine to munge the
file lines so that suitable values are extracted for comparison to $find.

    $tell = numeric( *BIG, $find, \&munge_line );
    $tell = alphabetic( *BIG, $find, \&munge_line );

A good example of this is the find_time() function. This is just an
implementation of the basic numeric algorithm similar to this.

    $tell = numeric ( *BIG, $epoch_seconds, \&get_epoch_seconds );

As the search is made the test lines are passed to the munging sub. This sub
needs to return a string or number that we can perform comparison on. In this
case the get_epoch_seconds sub looks for something that resembles a date
string, parses out the hh mm ss dd mm yyyy data, passes it to Time::Local for
conversion to epoch seconds and returns this number.

The optional arguments offered by Search::Dict to use dictionary order
(by removing non word/whitespace chars) and to ignore case are not directly
supported by File::SortedSeek, however they are easily implemented with a
simple munge function that does this:

    sub munge {
        local $_ = shift;
        s/[^\w\s]//g;  # removes non word/space chars for dictionary order
        return lc $_;  # makes comparison case insensistive
    }

=head2 get_epoch_seconds() - Convert a date string into epoch time

This function returns the epoch seconds represented by a date string in
the most common log file formats namely:

    Asctime format: Tue May 27 15:45:00 2008
    Apache format:  [21/May/2008:17:49:39 +1000]
    Squid format:   1012429341.115

The find_time() method uses this function internally. You can write your own
munge function if you need a different format. RTFS to see how but all the
function needs to do is take a data line and return a number that represents
the time. You could also use Date::Parse to do this for you.

=head2 find_time() - Seek to a specific time within the file

The find_time() function is an implementation of the basic numeric method as
discussed briefly above. You call it like:

    $tell = find_time( *LOG, 'Thu Jan  1 00:42:00 1970' );
    $tell = find_time( *LOG, $epoch_seconds );

You may use either a date string recognised by get_epoch_seconds() or just
epoch seconds directly.

=head2 get_between() - Getting lines from the middle of a file

Say you have a logfile and you want to get the log between one date and
another. You can simply use two calls to the find_time() to get the beginning
and end positions and then use get_between() to get the lines.

    # get the lines between two logfile dates
    $begin  = find_time( *LOG, $start );
    $end    = find_time( *LOG, $finish );
    # get lines as an array
    @lines = get_between( *LOG, $begin, $end );
    # get lines as an array reference
    $lines = get_between( *LOG, $begin, $end );

The get_between() method returns an array in list context as above and a
reference to an array in scalar context.

This function splits the lines based on a non system specific default record
separator heuristic. This is defined below:

    my $default_rec_sep = qw/\015\012|\015|\012/;

It should DWIM most of the time. You can override this on a per file basis
by passing the record separator to the get_between() function.

    @lines = get_between( *LOG, $begin, $end, $rec_sep );

Modifying $/ has no effect. Note that *the record separator is not returned*
in the array. As a result the returned array has effectively had every
element chomped.

Using the get_between() method you can efficiently get the lines at the
beginning of a file. Although you can just read in lines sequentially with
a while loop this requires that you test each line. If you can find the
end point using the find_time() numeric() or alphabetic() methods you
can the just get what you need. For large files many thousands of
unnecessary tests are avoided saving time. Using the example above
you simply set $begin to 0

    $begin  = 0;
    $end    = find_time( *LOG, $finish );
    @lines  = get_between( *LOG, $begin, $end );

You can similarly use get between to get all the lines from a specific point
up to the end of the file. The end is just the size of the file so:

    $begin = find_time( *LOG, $start );
    $end   = -s LOG;
    @lines = get_between( *LOG, $begin, $end );

=head2 get_last() - Get N lines from the end of a file.

This method does not depend on the file being sorted to work.
When you use the get_last() method the module estimates how many bytes at
the end of the file to read in. To make the estimate the module  multiplies
the default line length (80 chars) by the number of lines required and then
doubles it.

If it does not get sufficient lines on its first attempt it re-estimates
the line length from the actual data read in, re-calculates
the read, doubles it and then tries again. This algorithm is unlikely to
take more than 2 reads but if you have unusually long of short lines you may
get a small speed benefit by using the set_line_length() method to set the
average line length. The default is 80 chars per line. Setting the line length
close to the actual will also avoid reading a excessive quantity of data into
memory.

    # get last $n lines of any file as an array
    @lines = get_last( *BIG, $n )
    # or an array reference
    $lines = get_last( *BIG, $n )
    # change the input record separator from the default
    @lines = get_last( *BIG, $n, $rec_sep )

This function splits the lines based on a non system specific default record
separator heuristic. This is defined below:

    my $default_rec_sep = qw/\015\012|\015|\012/;

It should DWIM most of the time. You can override this on a per file basis
by passing the record separator to the get_between() function.
Example script in eg/tail.pl

=head1 EXPORT

Nothing by default. The following 5 methods are available for import:

    alphabetic()
    numeric()
    find_time()
    get_between()
    get_last()

You can import just the method(s) you want with:

    use File::SortedSeek qw(numeric);

or all 5  methods using the ':all' tag.

    use File::SortedSeek ':all';

=head1 OPTIONS and ACCESSORS

There are some options available via non exported function
calls. You will need to fully specify the name if you want to use these.

=head2 File::SortedSeek::error()

If a function returns undefined there has been an error. error() will
contain the text of the last error message or a null string if there
was no error.

=head2 File::SortedSeek::was_exact()

was_exact() will return true if an exact match was found. It will be
false if the match was in between or failed.

=head2 File::SortedSeek::set_cuddle()

set_cuddle() changes the default line returned for in between matches as
discussed above.

=head2 File::SortedSeek::set_no_cuddle()

set_no_cuddle() restores default behaviour

=head2 File::SortedSeek::set_descending()

By default ascending numerical order and alphabetical order are assumed.
This assumption can be reversed by calling set_descending() and reset
by calling set_ascending() We need to know the order to seek within the
file in the correct direction.

=head2 File::SortedSeek::set_ascending()

Reset sort order assumption to default.

=head2 File::SortedSeek::set_verbose()

Activate error messages. This is the default.

=head2 File::SortedSeek::set_silent()

Silence error messages

=head1 SPEED

Here is a table that demonstrates the advantage of using the binary search
algorithm.

    Num items    Lin av    Bin av   Lin:Bin
            2         1         1         1
            4         2         2         1
            8         4         3         1
           16         8         4         2
           32        16         5         3
           64        32         6         5
          128        64         7         9
          256       128         8        16
          512       256         9        28
         1024       512        10        51
         2048      1024        11        93
         4096      2048        12       170
         8192      4096        13       315
        16384      8192        14       585
        32768     16384        15      1092
        65536     32768        16      2048
       131072     65536        17      3855
       262144    131072        18      7281
       524288    262144        19     13797
      1048576    524288        20     26214

Even though there is an overhead involved with this search this is minor
as the number of tests required is so much less. Speed increases of 1000
of times are typical.

An OO interface slows things down by > 50% so is not used.

=head1 BUGS

All previously known ones have been removed.

=head1 CREDITS

Peter (Stig) Edwards for bugfixes and getting the refactoring started using
the cunning expedient of not just sending patches but also a refactored core.

Search::Dict for the basis of the code used to replace the original core
search function.

=head1 AUTHOR

(c) Dr James Freeman 2000-08 <airmedical [TA] gmail.com>
All rights reserved.

=head1 LICENSE

This package is free software and is provided "as is" without express or
implied warranty. It may be used, redistributed and/or modified under the
terms of the Artistic License 2.0. A copy is include in this distribution.

=head1 SEE ALSO

For details about the mystical significance of the number 42 and how it can
be applied to Life the Universe and everything see The Hitch Hiker's Guide
to the Galaxy 'trilogy' by the recently departed Douglas Adams.
