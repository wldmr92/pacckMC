classdef lookup_table
   
    properties 
        
        event = '';
        trans = [0, 0];
        rate = 0;
        psum = 0;
        
    end
    
    methods
        % constructor method
        function obj = lookup_table(event,transition)
            if nargin ~= 0
                m = size(event,1);
                n = size(event,2);
                obj(m,n) = obj;   
                for i=1:m
                   for j=1:n
                       obj(i,j).event = event{i,j};
                       obj(i,j).trans = transition(i,:);
                   end
                end
            end
        end
        
    end
    
end