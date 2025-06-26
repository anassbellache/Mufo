classdef test_simple_preprocess < matlab.unittest.TestCase
    methods (Test)
        function flatFieldCorrection(tc)
            proj = ones(4,5,6)*4;
            flat = ones(4,5)*5;
            dark = ones(4,5)*2;
            angles = linspace(0,pi,6);
            [corrected,sino,center] = mufo.simple_preprocess(proj,flat,dark,angles);
            tc.verifySize(corrected, size(proj));
            tc.verifySize(sino, [size(proj,2), size(proj,3), size(proj,1)]);
            tc.verifyClass(center,'double');
            expected = -log((4-2)/(5-2));
            tc.verifyEqual(corrected(1), expected, 'AbsTol',1e-12);
        end
        function reconstructSize(tc)
            proj = ones(4,5,6)*4;
            flat = ones(4,5)*5;
            dark = ones(4,5)*2;
            angles = linspace(0,pi,6);
            [~,sino,center] = mufo.simple_preprocess(proj,flat,dark,angles);
            volume = mufo.simple_reconstruct(sino, angles, center);
            tc.verifySize(volume, [size(sino,1), size(sino,1), size(sino,3)]);
        end
    end
end
