%% Clifford Algebra Geometric Library
% MSc. Edgar Macias Garcia (edgar.macias@cinvestav.mx)
% Dr. Julio Zamora Esquivel (julio.c.zamora.esquivel@intel.com )
% Prof. Eduardo José Bayro Corrochano (edb@gdl.cinvestav.mx)
% Centro de Investigación y Estudios Avanzados del Instituto Politécnico
% Nacional, Zapopan, México

%%

classdef blades

   properties
        basis = zeros(5,1);
        value = 0;
        grade = 0;
   end
    
   properties (Constant)
        ep = sym('ep', 'real');
        ec = sym('e', [5,1], 'real');
        e = [blades.ep; blades.ec];
        sb = [1; 1; 1; 1; -1];
        edf = sym('e', [5,5], 'real');
        edf3 = sym('e', [5,5,5], 'real');
        edf4 = sym('e', [5,5,5,5], 'real');
        edf5 = sym('e_1_2_3_4_5', 'real');
        ed = [blades.edf(1,2); blades.edf(1,3); blades.edf(1,4); blades.edf(1,5);
              blades.edf(2,3); blades.edf(2,4); blades.edf(2,5);
              blades.edf(3,4); blades.edf(3,5);
              blades.edf(4,5);
              blades.edf3(1,2,3); blades.edf3(1,2,4); blades.edf3(1,2,5); blades.edf3(1,3,4); blades.edf3(1,3,5); blades.edf3(1,4,5);
              blades.edf3(2,3,4); blades.edf3(2,3,5); blades.edf3(2,4,5);
              blades.edf3(3,4,5);
              blades.edf4(1,2,3,4); blades.edf4(1,2,3,5);
              blades.edf4(2,3,4,5);
              blades.edf5];
        base_ep = [0, 0, 0, 0, 0];
        base_e1 = [1, 0, 0, 0, 0];
        base_e2 = [0, 1, 0, 0, 0];
        base_e3 = [0, 0, 1, 0, 0];
        base_e4 = [0, 0, 0, 1, 0];
        base_e5 = [0, 0, 0, 0, 1];
        base = [blades.base_e1 + blades.base_e2;
                blades.base_e1 + blades.base_e3;
                blades.base_e1 + blades.base_e4;
                blades.base_e1 + blades.base_e5;
                blades.base_e2 + blades.base_e3;
                blades.base_e2 + blades.base_e4;
                blades.base_e2 + blades.base_e5;
                blades.base_e3 + blades.base_e4;
                blades.base_e3 + blades.base_e5;
                blades.base_e4 + blades.base_e5;
                blades.base_e1 + blades.base_e2 + blades.base_e3;
                blades.base_e1 + blades.base_e2 + blades.base_e4;
                blades.base_e1 + blades.base_e2 + blades.base_e5;
                blades.base_e1 + blades.base_e3 + blades.base_e4;
                blades.base_e1 + blades.base_e3 + blades.base_e5;
                blades.base_e1 + blades.base_e4 + blades.base_e5;
                blades.base_e2 + blades.base_e3 + blades.base_e4;
                blades.base_e2 + blades.base_e3 + blades.base_e5;
                blades.base_e2 + blades.base_e4 + blades.base_e5;
                blades.base_e3 + blades.base_e4 + blades.base_e5
                blades.base_e1 + blades.base_e2 + blades.base_e3 + blades.base_e4
                blades.base_e1 + blades.base_e2 + blades.base_e3 + blades.base_e5
                blades.base_e2 + blades.base_e3 + blades.base_e4 + blades.base_e5
                blades.base_e1 + blades.base_e2 + blades.base_e3 + blades.base_e4 + blades.base_e5];
            
        e0 = 0.5*(blades(5) - blades(4));
        ei = blades(4) + blades(5);
   end
   
   methods (Static)
       %Blades index function
       function b = basis_index(obj)
          b = [];
          for i = 1:size(obj.basis,1)
              if(obj.basis(i) == 1)
                  b(end+1) = i;
              end
          end
       end
       
       %Blades value
       function v = blade_param(obj, basis)                   
           v = 0;
           
           for i=1:size(obj,2)
              if(sum(obj(i).basis == basis) == 5)
                  v = obj(i).value;
                  break;
              end
           end
           
       end
       
       %Conformal point to vector
       function v = blade2vector(obj)
          
           v = [blades.blade_param(obj, blades(1).basis), blades.blade_param(obj, blades(2).basis), blades.blade_param(obj, blades(3).basis)];
           
       end
       
       %Extract one dimensionals from multivector
       function v = onedim(obj)
          
           v = blades.blade_param(obj, blades(0).basis)*blades(0) + blades.blade_param(obj, blades(1).basis)*blades(1) + blades.blade_param(obj, blades(2).basis)*blades(2) + blades.blade_param(obj, blades(3).basis)*blades(3) + blades.blade_param(obj, blades(4).basis)*blades(4) + blades.blade_param(obj, blades(5).basis)*blades(5); 
       end
       
       %Operations with conformal points
       function p = confpsub(obj1, obj2)
          
           %Conformals to vectors
           p1 = blades.blade2vector(obj1);
           p2 = blades.blade2vector(obj2);
           
           %Do substraction
           p = blades.euc2confpoint(p1 - p2);
           
       end
       
       function p = confpsum(obj1, obj2)
          
           %Conformals to vectors
           p1 = blades.blade2vector(obj1);
           p2 = blades.blade2vector(obj2);
           
           %Do substraction
           p = blades.euc2confpoint(p1 + p2);
           
       end
       
       function n = confpnorm(obj1)
          
           %Conformals to vectors
           p = blades.blade2vector(obj1);
           
           %Do substraction
           n = sqrt(p*p');
           
       end
       
       %Sign calculation according position of basis 
       function s = dot_sign(obj1, obj2)
          
           %Get index of blades
           ba = blades.basis_index(obj1);
           bb = blades.basis_index(obj2);
           s = 1;
           
           %Check basis multiplication
           for i = 0:size(ba,2)-1
               d = 0;
               %Existing basis
               e = ba(end - i);
               %Checkout position
               for j = 1:size(bb,2)
                   if(bb(j) < e)
                       d = d + 1;
                   end
               end
                   
               %Check number of positions
               if(mod(d,2) == 1)
                   s = s*(-1);
               end
           end
       end
       
       %function d = confdistance(a , b)
          %Calculate distance between points
          
           
       %end
       
       %Blade simplification function
       function obj3 = simplify_blade(obj)
           
          if(size(obj,2) == 1)
              obj3 = obj;
          else
              obj2 = obj(1);
              obj(1) = [];
              obj3 = [];

              %Check common basis
              while(size(obj,2) > 0)
                  oba = [];
                  for j=1:size(obj,2)
                      if(sum(obj2(end).basis == obj(j).basis)==5)
                         obj2(end).value = obj2(end).value + obj(j).value;
                         oba(end+1) = j;
                      end
                  end

                  %Clear basis
                  for j=0:size(oba,2)-1
                     obj(oba(end - j)) = []; 
                  end

                  %Append new axis
                  if(size(obj,2) > 0)
                      obj2(end+1) = obj(1);
                      obj(1) = [];
                  end
              end

              %Erase basis with zero value
              for i = 1:size(obj2,2)
                 if(abs(obj2(i).value) > 0.0001)
                     if(size(obj3) == 0)
                         obj3 = obj2(i);
                     else
                         obj3(end + 1) = obj2(i);
                     end
                 end
              end
          end
       end
       
       %Euclidean to conformal conversion
       function obj = euc2confpoint(v)      
           %Check size
           if(size(v,2) == 3)
              obj = v(1)*blades(1) + v(2)*blades(2) + v(3)*blades(3) + 0.5*(v*v')*(blades(4) + blades(5)) + 0.5*(blades(5)-blades(4));          
           else
               disp('Invalid input')
               obj = 0;
           end        
       end
       
       %Euclidean vector to conformal point
       function v = eucp2confp(v_c)
          %Point to vector form
          v_e = blades.blade2vector(v_c);
          
          %Vector to conformal point
          v = blades.euc2confpoint(v_e);
           
       end
       
       %Conformal point to Euclidean point
       function obj = confpoint2euc(v)  
              obj = blades.blade_param(v, blades(1).basis)*blades(1) + blades.blade_param(v, blades(2).basis)*blades(2) + blades.blade_param(v, blades(3).basis)*blades(3);          
       end
       
       function obj = euc2conf(v)      
           %Check size
           if(size(v,2) == 3)
              obj = v(1)*blades(1) + v(2)*blades(2) + v(3)*blades(3);          
           else
               disp('Invalid input')
               obj = 0;
           end        
       end
       
       %Square norm
       function n = norm(obj)
           %Get square norm
           N = blades.simplify_blade(obj*obj);
           n = sqrt(abs(N.value));
       end
       
       function n = eucnorm(obj)
          
           %Get componentes
           x = blades.blade_param(obj, blades(1).basis);
           y = blades.blade_param(obj, blades(2).basis);
           z = blades.blade_param(obj, blades(3).basis);
           
           %Get squared norm
           n = sqrt(x*x + y*y + z*z);
       end
   end
    
   methods
       
        %Blades construtor
        function obj = blades(x)
            if(size(x,2) == 1)
                if(x==0)
                    obj.grade = 0;
                    obj.value = 1;  
                else
                    obj.basis(x) = 1;
                    obj.grade = 1;
                    obj.value = 1;
                end
            elseif (size(x,2) == 2)
                a = zeros(5,1); a(x(1)) = 1;
                b = zeros(5,1); b(x(2)) = 1;
                obj.basis = a | b;
                obj.grade = 2;
                if(x(2)>x(1))
                    obj.value = 1;
                else
                    obj.value = -1;
                end
            elseif (size(x,2) == 3)
                a = zeros(5,1); a(x(1)) = 1;
                b = zeros(5,1); b(x(2)) = 1;
                c = zeros(5,1); c(x(3)) = 1;
                obj.basis = a | b | c;
                obj.grade = 3;
                obj.value = 1;
            elseif (size(x,2) == 5)
                obj.basis = x';
                obj.grade = sum(x);
                obj.value = 1;
            end
        end
        
        %Sum operation
        function obj = plus(obj1, obj2)
            %Sum by a zero constant
            if(isa(obj1, 'double'))
                obj_a = blades(0);
                if(size(obj1,1) == 0 && size(obj1,2) == 0)
                    obj_a.value = 0;
                else
                    obj_a.value = obj1;
                end
                obj = obj2;
                obj1 = obj2;
                obj2 = obj_a;
            elseif(isa(obj2, 'double'))
                obj_a = blades(0);
                if(size(obj2,1) == 0 && size(obj2,2) == 0)
                    obj_a.value = 0;
                else
                    obj_a.value = obj2;
                end
                obj = obj1;
                obj2 = obj_a;
            else
                obj = obj1;
            end
            
            for i = 1:size(obj2,2)
                cond = true;
                for j = 1:size(obj1,2)
                    %Check basis
                    if(sum(obj1(j).basis == obj2(i).basis)==5)
                        obj(j).value = obj(j).value + obj2(i).value;
                        cond = false;
                    end
                end
                
                %Add new basis
                if(cond)
                   obj(end+1) = obj2(i); 
                end
            end
            
        end
        
        %Rest operation
        function obj = minus(obj1, obj2)
            obj = obj1;
 
            for i = 1:size(obj2,2)
                cond = true;
                for j = 1:size(obj1,2)
                    %Check basis
                    if(sum(obj1(j).basis == obj2(i).basis)==5)
                        obj(j).value = obj(j).value - obj2(i).value;
                        cond = false;
                    end
                end
                
                %Add new basis
                if(cond)
                   obj(end+1) = -1*obj2(i); 
                end
            end
            
        end
        
        %Clifford product operation
        function obj = mtimes(a,b)
            if(isa(a, 'double') && isa(b, 'blades'))
                obj = b;
                for i=1:size(b,2)
                    obj(i).value = b(i).value*a;
                end
            elseif(isa(b, 'double') && isa(a, 'blades'))
                obj = a;
                for i=1:size(a,2)
                    obj(i).value = a(i).value*b;
                end
            elseif(isa(a, 'blades') && isa(b, 'blades'))
                for i=1:size(a,2)
                    for j=1:size(b,2)
                        v = xor(a(i).basis, b(j).basis)';
                        s = and(a(i).basis, b(j).basis)';
                        m = blades.dot_sign(a(i), b(j));
                        
                        for k=1:size(s,2)
                           if(s(k) ~= 0)
                              m = m*blades.sb(k);
                           end
                        end
                        if(i == 1 && j ==1)
                            obj2 = blades(v');
                            obj2.value = a(i).value*b(j).value*m;
                            obj2.grade = sum(obj2.basis);
                        else
                            obj2(end+1) = blades(v');
                            obj2(end).value = a(i).value*b(j).value*m;
                            obj2(end).grade = sum(obj2(end).basis);
                        end
                    end
                end
                obj = obj2;
            elseif(isa(a, 'sym') && isa(b, 'blades'))
                obj = b;
                for i=1:size(b,2)
                    obj(i).value = b(i).value*a;
                end
            elseif(isa(b, 'sym') && isa(a, 'blades'))
                obj = a;
                for i=1:size(b,2)
                    obj(i).value = a(i).value*b;
                end    
            else
                disp('Invalid operation')
            end
            
            %Simplify array
            obj = blades.simplify_blade(obj);
        end
        
        %Wedge product operation
        function obj = mpower(obj1, obj2)
            obj = 0.5*(obj1*obj2 - obj2*obj1);
            obj = blades.simplify_blade(obj);
            
        end
        
        %Dot product operation
        function obj = times(obj1, obj2)
            obj = 0.5*(obj1*obj2 + obj2*obj1);
            obj = blades.simplify_blade(obj);
        end
        
        %Blades display function
        function disp(obj)
            %One component
            if(size(obj,2) == 1)
                r = 0;
                if(obj.grade == 1 || obj.grade == 0)
                    if(sum(obj.basis) == 0)
                        r = r + obj.value*blades.e(1);
                    else
                        for i = 1:5
                            if(obj.basis(i) == 1)
                               r = r + obj.value*blades.e(i+1);
                            end
                        end
                    end
                else
                    for i = 1:size(blades.base,1)
                        if(sum(obj.basis' == blades.base(i,:)) == 5)
                            r = r + obj.value*blades.ed(i);
                        end 
                    end
                end
                disp(vpa(r,6));
            %Multivector
            else
                r = 0;
                for i = 1:size(obj,2)
                    if(obj(i).grade < 2)
                        if(sum(obj(i).basis) == 0)
                            r = r + obj(i).value*blades.e(1);
                        else
                            for j = 1:5
                                if(obj(i).basis(j) == 1)
                                   r = r + obj(i).value*blades.e(j+1);
                                end
                            end
                        end
                    else
                        for k = 1:size(blades.base,1)
                            if(sum(obj(i).basis' == blades.base(k,:)) == 5)
                                r = r + obj(i).value*blades.ed(k);
                            end 
                        end
                    end
                end
                disp(vpa(r,6));
            end
        end
    end
end