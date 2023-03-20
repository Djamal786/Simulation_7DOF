robot = importrobot('iiwa14.urdf');
test = struct('JointPositions',cell(1,7))

qPositions = num2cell(qSim7);
for i = 1:200
     config = homeConfiguration(robot)  
     for j = 1:7
        config(j).JointPosition = qSim7(j,i)
     end
    show(robot,config);
    drawnow
end
  