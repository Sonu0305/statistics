## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{x} =} probit (@var{p})
##
## Probit transformation
##
## Return the probit (the quantile of the standard normal distribution) for
## each element of @var{p}.
##
## @seealso{logit}
## @end deftypefn

function x = probit (p)

  if (nargin != 1)
    print_usage ();
  endif

  if (! isnumeric (p) && ! islogical (p))
    error ("probit: P must be a numeric or logical array.");
  endif

  if (! isfloat (p))
    p = double (p);
  endif

  x = -sqrt (2) * erfcinv (2 * p);

endfunction

## Test output
%!assert_equal (probit ([-1, 0, 0.5, 1, 2]), [NaN, -Inf, 0, Inf, NaN])
%!assert_equal (probit ([0.2, 0.99]), norminv ([0.2, 0.99]))

## Test input validation
%!error probit ()
%!error probit (1, 2)
%!error <probit: P must be a numeric or logical array.> probit ("a")
%!error <probit: P must be a numeric or logical array.> probit ({0.5})

## Test logical and integer casting
%!test
%! assert_equal (probit (true), Inf);
%! assert_equal (probit (false), -Inf);
%! assert_equal (probit (int8(0)), -Inf);

## Test on large 3D numeric array (huge complex inputs)
%!test
%! p = ones (50, 50, 50) * 0.5;
%! x = probit (p);
%! assert_equal (size (x), [50, 50, 50]);
%! assert_equal (x, zeros (50, 50, 50));
