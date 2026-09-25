function paths = cfar_setup()
%CFAR_SETUP  Put every part of this project on the MATLAB path and report
%   where the datasets live. Run once per session (it is idempotent).
%
%   paths = CFAR_SETUP()
%
%   Returns a struct of absolute paths:
%     .root, .common, .comparison, .results, .ssdd, .ssddJpeg, .ssddXml,
%     .mstar, .hrsid, .hrsidImages, .hrsidAnnFile, and one field per
%     detector folder.
%
%   Layout this enforces:
%     _common/          shared front end, solvers, metrics  (all detectors)
%     _datasets/SSDD/   1160 images + VOC annotations, self-contained copy
%     _datasets/MSTAR/  12 unique MSTAR chips, deduplicated
%     _comparison/      the Phase 2 cross-detector sweep and its outputs
%     CFAR <Dist>/      one folder per detector: Params, TLog, Floating,
%                       demo drivers, and legacy/ holding the originals
%
%   Detector code lives in its own folder; only genuinely shared
%   infrastructure lives in _common.

root = fileparts(mfilename('fullpath'));

paths = struct();
paths.root       = root;
paths.common     = fullfile(root, '_common');
paths.comparison = fullfile(root, '_comparison');
paths.results    = fullfile(root, '_comparison', 'Results');
paths.figures    = fullfile(root, '_comparison', 'Figures');

paths.weibull    = fullfile(root, 'CFAR_Weibull');
paths.lognormal  = fullfile(root, 'CFAR Lognormal');
paths.gengamma   = fullfile(root, 'CFAR Generalized Gamma');
paths.g0         = fullfile(root, 'CFAR_G0');
paths.burr       = fullfile(root, 'CFAR_Burr');
paths.k          = fullfile(root, 'CFAR K');

paths.ssdd       = fullfile(root, '_datasets', 'SSDD');
paths.ssddJpeg   = fullfile(paths.ssdd, 'JPEGImages');
paths.ssddXml    = fullfile(paths.ssdd, 'Annotations');
paths.mstar      = fullfile(root, '_datasets', 'MSTAR');
paths.sarfish    = fullfile(root, 'sarfish_sample', 'crops');

paths.hrsid          = fullfile(root, 'HRSID');
paths.hrsidImages    = fullfile(paths.hrsid, 'images');
paths.hrsidAnnFile   = fullfile(paths.hrsid, 'annotations', 'train_test2017.json');

addpath(paths.common, paths.comparison, ...
        paths.weibull, paths.lognormal, paths.gengamma, paths.g0, paths.burr, paths.k);

for f = {'results', 'figures'}
    if ~isfolder(paths.(f{1}))
        mkdir(paths.(f{1}));
    end
end

if ~isfolder(paths.ssddJpeg)
    warning('cfar_setup:NoSSDD', ...
        ['SSDD images not found at %s. The Phase 2 comparison needs them ' ...
         '(they carry the ground-truth boxes Pd/Pfa are measured against).'], ...
        paths.ssddJpeg);
end

if ~isfolder(paths.hrsidImages)
    warning('cfar_setup:NoHRSID', ...
        ['HRSID images not found at %s. The HRSID comparison sweep ' ...
         '(run_comparison_hrsid) needs them.'], paths.hrsidImages);
end
end
