classdef testUmv < matlab.unittest.TestCase
    
    
    methods(Test)
        
        function testEstimatedInput(testCase)
            % Verify that alpha rad is less than 0.1rad after 10s
            %endPoint = testCase.AirframeBusData.alpha_rad.Data(end);
            %testCase.verifyLessThan(endPoint, 0.1);
        end
        
        function testEstimatedState(testCase)
            
        end
        
    end
end