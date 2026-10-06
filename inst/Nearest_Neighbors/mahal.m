## Copyright (C) 2015 Lachlan Andrew
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
## @deftypefn  {statistics} {@var{d} =} mahal (@var{y}, @var{x})
##
## Mahalanobis' D-square distance.
##
## Return the Mahalanobis' D-square distance of the points in
## @var{y} from the distribution implied by points @var{x}.
##
## Specifically, it uses a Cholesky decomposition to set
##
## @example
##  answer(i) = (@var{y}(i,:) - mean (@var{x})) * inv (A) *
##              (@var{y}(i,:)-mean (@var{x}))'
## @end example
##
## where A is the covariance of @var{x}.
##
## The data @var{x} and @var{y} must have the same number of components
## (columns), but may have a different number of observations (rows).
##
## @end deftypefn

function retval = mahal (y, x)

  if (nargin != 2)
    print_usage ();
  endif

  if (! isfloat (x))
    x = double (x);
  endif
  if (! isfloat (y))
    y = double (y);
  endif

  if (! isreal (x) || ! isreal (y))
    error ("mahal: X and Y must not be complex arrays.");
  endif

  if (! ismatrix (x) || ! ismatrix (y) || size (x, 2) != size (y, 2))
    error ("mahal: X and Y must have an identical number of columns.");
  endif

  [xr, xc] = size (x);
  [yr, yc] = size (y);

  xm = mean (x, 1);

  ## Center data by subtracting mean of x
  x = bsxfun (@minus, x, xm);
  y = bsxfun (@minus, y, xm);

  w = (x' * x) / (xr - 1);

  retval = sumsq (y / chol (w), 2);

endfunction


## Test input validation
%!error mahal ()
%!error mahal (1, 2, 3)
%!error <mahal: X and Y must not be complex arrays.> mahal ([1+1i], [2+2i])
%!error <must have an identical number of columns> mahal (ones (2, 2), ones (2, 2, 2))
%!error <must have an identical number of columns> mahal (ones (2, 2), ones (2, 3))

%!test
%! X = [1 0; 0 1; 1 1; 0 0];
%! assert_equal (mahal (X, X), [1.5; 1.5; 1.5; 1.5], 10*eps)
%! assert_equal (mahal (X, X+1), [7.5; 7.5; 1.5; 13.5], 10*eps)

%!assert_equal (mahal ([true; true], [false; true]), [0.5; 0.5], eps)
