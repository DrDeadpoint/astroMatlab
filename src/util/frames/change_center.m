function traj_out = change_center(traj_in, new_center)
    currentFrame = traj_in.system_model.frame;
    if ~strcmp(currentFrame(end-2:end),'rot')
        error('center change currently only supported in rotating frame')
    end
    mu = traj_in.system_model.char.mu;
    B1 = 0;
    P1 = -mu;
    P2 = 1-mu;
    current_cnter = currentFrame(1:2);
    switch current_cnter
        case 'B1'
            cx = B1;
        case 'P1'
            cx = P1;
        case 'P2'
            cx = P2;
    end
    switch new_center
        case 'B1'
            nx = B1;
        case 'P1'
            nx = P1;
        case 'P2'
            nx = P2;
    end
    diffx = nx - cx;
    traj_out = traj_in;
    traj_out.pos = traj_out.pos.change_unit('nd_l', traj_out.system_model);
    traj_out.pos.value(1,:) = traj_out.pos.value(1,:) - diffx;
    traj_out.system_model.frame = [new_center currentFrame(3:end)];
end