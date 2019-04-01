function [n_alpha, n_single, p_alpha, alphas, separable_fraction]...
    = SeparabilityAnalysis(X, varargin)
%Performs standard analysis of separability and produces standard plots. 
%
%Inputs:
%   X  - is a data matrix with one data point in each row.
%   Optional arguments in varargin form Name, Value pairs. Possible names:
%       'ConditionalNumber' - a positive real value used to select the top
%           princinpal components. We consider only PCs with eigen values
%           which are not less than the maximal eigenvalue divided by
%           ConditionalNumber Default value is 10.
%       'ProjectOnSphere' - a boolean value indicating if projecting on a
%           sphere should be performed. Default value is true.
%       'Alphas' - a real vector, with alpha range, the values must be
%           within (0,1) interval. Default is [0.6,0.62,...,0.98].
%       'ProducePlots' - a boolean value indicating if the standard plots
%           needs to be drawn. Default is true.
%
%Outputs:
%   n_alpha - effective dimension profile as a function of alpha
%   n_single - a single estimate for the effective dimension 
%   p_alpha - distributions as a function of alpha, matrix with columns
%       corresponding to the alpha values, and with rows corresponding to
%       objects. 
%   separable_fraction - separable fraction of data points as a function of
%       alpha
%   alphas - alpha values
% 

    % Define default values of optional parameters
    ConditionalNumber = 10;
    ProjectOnSphere = 1;
    alphas = 0.6:0.02:0.98;
    ProducePlots = 1;
    % Read optional parameters
    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'ConditionalNumber')
            ConditionalNumber = varargin{i + 1};
        end
        if strcmpi(varargin{i},'ProjectOnSphere')
            ProjectOnSphere = varargin{i + 1};
        end
        if strcmpi(varargin{i},'Alphas')
            alphas = varargin{i + 1};
        end
        if strcmpi(varargin{i},'ProducePlots')
            ProducePlots = varargin{i + 1};
        end
    end
    % Sanity check of parameters
    if ~isnumeric(ConditionalNumber) || ConditionalNumber < 0
        error('ConditionalNumber must be a positive real number');
    end
    % Transform expected boolean values into Boolean
    if ~islogical(ProjectOnSphere)
        if isnumeric(ProjectOnSphere)
            ProjectOnSphere = ProjectOnSphere ~= 0;
        else
            error('ProjectOnSphere must be Boolean value "true" or number 1');
        end
    end
    if ~islogical(ProducePlots)
        if isnumeric(ProducePlots)
            ProducePlots = ProducePlots ~= 0;
        else
            error('ProducePlots must be Boolean value "true" or number 1');
        end
    end
    alphas = alphas(:)';
    if ~isnumeric(alphas) || sum(alphas <= 0) > 0 || sum(alphas >= 1) > 0
        error(['"Alphas" must be a real vector, with alpha range,'...
            ' the values must be within (0,1) interval']);
    end
    
    % Sort alphas
    alphas = sort(alphas);
    
    % Preprocess data
    Xp = preprocessing(X, 1, 1, 1, ProjectOnSphere,...
        'ConditionalNumber',ConditionalNumber);
    % Check separability
    [separable_fraction, p_alpha] = checkSeparabilityMultipleAlpha(Xp, alphas);
    
    % Calculate mean of fraction of separable for each alpha.
    py_mean = mean(p_alpha,2);

    [n_alpha, n_single] = dimension_uniform_sphere(py_mean, alphas);

    alpha_ind_selected = find(n_single==n_alpha);

    % Draw figures if requested
    if ProducePlots
        % Define the minimal and maximal dimensions for theoretical graph with
        % two dimensions in each side.
        n_min = floor(min(n_alpha)) - 2;
        n_max = floor(max(n_alpha)+0.8) + 2;
        %n_min must be at least 1
        if n_min < 1
            n_min = 1;
        end
        ns = n_min:n_max;
        % Empirical graph of dimension vs alpha
        subplot(1,3,1);
        plot(alphas, n_alpha, 'ko-');
        hold on;
        plot(alphas(alpha_ind_selected), n_single, 'rx', 'MarkerSize', 20);
        xlabel('Alpha', 'FontSize', 16);
        ylabel('Effective dimension', 'FontSize', 16);
        
        % Distribution of fraction of inseparable points for reference alpha
        % Define number of bins
        nbins = floor(size(X,1)/200);
        if nbins<20
            nbins = 20;
        end
        subplot(1,3,2);
        hist(p_alpha(alpha_ind_selected,:),nbins);
        xlabel(sprintf('inseparability prob for alpha=%2.2f',...
            alphas(alpha_ind_selected)), 'FontSize', 16);
        ylabel('Number of values', 'FontSize', 16);
        
        % Theoretical curves
        subplot(1,3,3);
        semilogy(alphas, py_mean, 'bo-', 'LineWidth', 3);
        hold on;
        xlabel('Alpha', 'FontSize', 16);
        ylabel('Mean fraction of inseparabile', 'FontSize', 16);
        title(sprintf('Theor.curves for n=%i:%i',n_min,n_max), 'FontSize', 16);
        % Calculate and draw theoretical curves
        pteor = probability_unseparable_sphere(alphas,ns(:));
        semilogy(alphas, pteor, '-');
        set(gcf,'Position',[105, 161, 1115, 344]);
    end
end
