


function fxns = plot_fxns()

    fxns.get_empty_vector = @get_empty_vector;


end

function plot_struct = get_plot_structure()
    plot_struct.name.title = '';
    plot_struct.name.y = '';
    plot_struct.name.x = '';
    plot_struct.length = NaN;
    plot_struct.color = '';
end

function vector = get_empty_vector(length)
    vector = zeros(length, 1);

end