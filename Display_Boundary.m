function Display_Boundary(XL, XR, YL, YR)
    Linevalue=1;
    plot([XL XL], [YL YR], '-k', 'LineWidth',Linevalue);
    plot([XR XR], [YL YR], '-k', 'LineWidth',Linevalue);
    plot([XL XR], [YL YL], '-k', 'LineWidth',Linevalue);
    plot([XL XR], [YR YR], '-k', 'LineWidth',Linevalue);
end