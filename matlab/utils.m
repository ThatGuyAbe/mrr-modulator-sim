function U = utils()
% small helpers shared by the step scripts, as function handles:
%   U.tee(fid, fmt, ...)     print to the console and to an open summary file
%   U.yesno(b)               'yes' or 'NO'
%   U.out_folder(S, name)    output folder <S.out_root>/<name>, created if needed
%   U.save_png(fig, fname)   save a quick-look figure and close it
%   U.cross_at(x, y, lev)    first crossing of level lev along (x, y)
%   U.load_targets(S)        dw_opt, dw_lin and the sweeps from step 1
%   U.start_pool()           open a parallel pool, returns the number of workers
    U.tee          = @tee;
    U.yesno        = @yesno;
    U.out_folder   = @out_folder;
    U.save_png     = @save_png;
    U.cross_at     = @cross_at;
    U.load_targets = @load_targets;
    U.start_pool   = @start_pool;
end

function tee(fid, varargin)
    fprintf(1, varargin{:});
    if fid > 0, fprintf(fid, varargin{:}); end
end

function s = yesno(b)
    if b, s = 'yes'; else, s = 'NO'; end
end

function d = out_folder(S, name)
    d = fullfile(S.out_root, name);
    if ~exist(d, 'dir'), mkdir(d); end
end

function save_png(fig, fname)
    try
        print(fig, fname, '-dpng', '-r150');
    catch err
        fprintf('quick-look figure not saved (%s)\n', err.message);
    end
    close(fig);
end

function x0 = cross_at(x, y, lev)
% linear interpolation between the two samples on either side of the crossing
    s = sign(y(:) - lev);
    k = find(s(1:end-1) ~= s(2:end), 1);
    if isempty(k), x0 = NaN; return; end
    x0 = x(k) + (lev - y(k)) * (x(k+1) - x(k)) / (y(k+1) - y(k));
end

function T = load_targets(S)
% T.dw_opt(iQ, iP, iR) for S.Q_values x S.powers x S.datarates (GHz), T.dw_veo the
% veo maximum on the same curves, T.dw_lin(iQ), and T.get(Q, P_mW, rate) for one value
    f = fullfile(S.out_root, 'step1_sweeps', 'results.mat');
    if exist(f, 'file') ~= 2
        error('run step1_sweeps.m first (%s not found)', f);
    end
    L = load(f);
    T.Q_values = L.Q_values;  T.powers = L.powers;  T.datarates = L.datarates;
    T.dw_opt = L.dw_opt;  T.dw_veo = L.dw_veo;  T.dw_lin = L.dw_lin;
    T.df_grid = L.df_grid;  T.oma_sw = L.oma_sw;  T.veo_sw = L.veo_sw;
    T.get = @(Q, Pm, r) L.dw_opt(find(L.Q_values == Q, 1), find(abs(L.powers - Pm) < 1e-9, 1), ...
                                 find(abs(L.datarates - r) < 1e-9, 1));
    T.get_veo = @(Q, Pm, r) L.dw_veo(find(L.Q_values == Q, 1), find(abs(L.powers - Pm) < 1e-9, 1), ...
                                     find(abs(L.datarates - r) < 1e-9, 1));
end

function nW = start_pool()
% tries an open pool, then a process pool ('Processes', called 'local' before
% r2022b), then a thread pool. if none starts, parfor runs serially with the same
% results. in octave parfor is a plain loop.
    if exist('OCTAVE_VERSION', 'builtin') ~= 0
        nW = 1;
        return;
    end
    pool = gcp('nocreate');
    if isempty(pool), pool = try_pool('Processes'); end
    if isempty(pool), pool = try_pool('local'); end
    if isempty(pool), pool = try_pool('Threads'); end
    if isempty(pool)
        try
            ps = parallel.Settings;
            ps.Pool.AutoCreate = false;        % so parfor does not try again
        catch
        end
        warning(['No parallel pool could be started. parfor loops run serially ' ...
                 '(same results, about N times slower for N cores).']);
        nW = 1;
        return;
    end
    files = {'sim_mrr_eye.m', 'sim_mrr_fulltrace.m', 'sim_mrr_cw_sweep.m', 'mod_rk4_step.m', ...
             'run_case.m', 'ring_params.m', 'make_drive.m'};
    try
        addAttachedFiles(pool, files);     % process pools only, thread pools share the session
    catch
    end
    nW = pool.NumWorkers;
    fprintf('parallel pool: %s, %d workers\n', class(pool), nW);
end

function pool = try_pool(profile)
    pool = [];
    try
        pool = parpool(profile);
    catch err
        fprintf('parpool(''%s'') failed: %s\n', profile, err.message);
    end
end
