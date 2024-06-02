function [densities, varargout] = solver(S, bc, rhs, eps, varargin)
%LAP3D.SOLVER: solve Laplace boundary value problems using integral 
% equations
%
% Syntax for Dirichlet problems
%   [densities] = lap3d.solver(S, 'dir', rhs, eps);
%   [densities] = lap3d.solver(S, 'dir', rhs, eps, rep_pars);
%   [densities] = lap3d.solver(S, 'dir', rhs, eps, rep_pars, opts);
%
% Syntax for Neumann problems
%   [densities] = lap3d.solver(S, 'neu', rhs, eps);
%   [densities] = lap3d.solver(S, 'neu', rhs, eps, opts);
%
    switch lower(bc)
      case {'dir', 'd', 'dirichlet'}
        if nargin < 5
          rep_params = [1; 1];
          opts = [];
        elseif nargin < 6
          if isa(varargin{1}, 'double')
            rep_params = varargin{1};
          else
            error('LAP3D.SOLVER: rep_params, should be a double array');
          end
          opts = [];
        else 
          if isa(varargin{1}, 'double')
            rep_params = varargin{1};
          else
            error('LAP3D.SOLVER: rep_params, should be a double array');
          end
          opts = varargin{2};
        end
        [densities, varargout{2:nargout}] = lap3d.dirichlet.solver(S, rhs, eps, ... 
                                    rep_params, opts);

      case {'neu', 'n', 'neumann'}
        if nargin < 5
          opts = [];
        else
          opts = varargin{1};
        end
        [densities, varargout{2:nargout}] = lap3d.neumann.solver(S, rhs, eps, opts);
    end
    
end