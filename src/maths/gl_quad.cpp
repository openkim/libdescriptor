#include "maths/gl_quad.hpp"
#include <vector>

// gauss-Lebedev quadrature for size 100, same as taken by dscribe
// TODO: higher accuracy quadrature?

std::vector<double> get_gl_weights() {
    return {7.34634490505672E-4, 0.001709392653518105, 0.002683925371553482, 0.003655961201326375,
            0.00462445006342212, 0.00558842800386552, 0.00654694845084532, 0.007499073255464712,
            0.008443871469668971, 0.00938041965369446, 0.01030780257486897, 0.01122511402318598,
            0.012131457662979497, 0.0130259478929715423, 0.0139077107037187727, 0.014775884527441302,
            0.015629621077546003, 0.01646808617614521, 0.017290460568323582, 0.018095940722128117,
            0.018883739613374905, 0.019653087494435306, 0.020403232646209433, 0.021133442112527642,
            0.02184300241624739, 0.022531220256336273, 0.023197423185254122, 0.023840960265968206,
            0.02446120270795705, 0.02505754448157959, 0.02562940291020812, 0.02617621923954568,
            0.026697459183570963, 0.02719261344657688, 0.027661198220792388, 0.028102755659101173,
            0.028516854322395098, 0.028903089601125203, 0.029261084110638277, 0.02959048805991264,
            0.029890979593332831, 0.03016226510516914, 0.03040407952645482, 0.03061618658398045,
            0.03079837903115259, 0.030950478850490988, 0.031072337427566517, 0.031163835696209907,
            0.03122488425484936, 0.03125542345386336, 0.031255423453863357, 0.03122488425484936,
            0.0311638356962099068, 0.031072337427566517, 0.030950478850490988, 0.03079837903115259,
            0.030616186583980449, 0.03040407952645482, 0.0301622651051691449, 0.02989097959333283,
            0.029590488059912643, 0.029261084110638277, 0.028903089601125203, 0.0285168543223951,
            0.02810275565910117, 0.02766119822079239, 0.02719261344657688, 0.02669745918357096,
            0.02617621923954568, 0.025629402910208116, 0.02505754448157959, 0.024461202707957053,
            0.02384096026596821, 0.023197423185254122, 0.0225312202563362727, 0.021843002416247386,
            0.02113344211252764, 0.020403232646209433, 0.019653087494435306, 0.0188837396133749046,
            0.018095940722128117, 0.017290460568323582, 0.016468086176145213, 0.015629621077546003,
            0.0147758845274413, 0.013907710703718773, 0.01302594789297154, 0.0121314576629795,
            0.01122511402318598, 0.01030780257486897, 0.00938041965369446, 0.008443871469668971,
            0.007499073255464712, 0.00654694845084532, 0.005588428003865515, 0.00462445006342212,
            0.003655961201326375, 0.002683925371553482, 0.001709392653518105, 7.3463449050567E-4};
}

std::vector<double> get_orig_gl_grid(){
    return {-0.999713726773441234, -0.998491950639595818, -0.996295134733125149, -0.99312493703744346,
            -0.98898439524299175, -0.98387754070605702, -0.97780935848691829, -0.97078577576370633,
            -0.962813654255815527, -0.95390078292549174, -0.94405587013625598, -0.933288535043079546,
            -0.921609298145333953, -0.90902957098252969, -0.895561644970726987, -0.881218679385018416,
            -0.86601468849716462, -0.849964527879591284, -0.833083879888400824, -0.815389238339176254,
            -0.79689789239031448, -0.77762790964949548, -0.757598118519707176, -0.736828089802020706,
            -0.715338117573056447, -0.69314919935580197, -0.670283015603141016, -0.64676190851412928,
            -0.622608860203707772, -0.59784747024717872, -0.57250193262138119, -0.546597012065094168,
            -0.520158019881763057, -0.493210789208190934, -0.465781649773358042, -0.437897402172031513,
            -0.409585291678301543, -0.380872981624629957, -0.351788526372421721, -0.322360343900529152,
            -0.292617188038471965, -0.26258812037150348, -0.23230248184497397, -0.201789864095735997,
            -0.171080080538603275, -0.140203137236113973, -0.109189203580061115, -0.0780685828134366367,
            -0.046871682421591632, -0.015628984421543083, 0.0156289844215430829, 0.046871682421591632,
            0.078068582813436637, 0.109189203580061115, 0.140203137236113973, 0.171080080538603275,
            0.201789864095735997, 0.23230248184497397, 0.262588120371503479, 0.292617188038471965,
            0.322360343900529152, 0.351788526372421721, 0.380872981624629957, 0.409585291678301543,
            0.437897402172031513, 0.465781649773358042, 0.49321078920819093, 0.520158019881763057,
            0.546597012065094168, 0.572501932621381191, 0.59784747024717872, 0.622608860203707772,
            0.64676190851412928, 0.670283015603141016, 0.693149199355801966, 0.715338117573056447,
            0.736828089802020706, 0.75759811851970718, 0.77762790964949548, 0.79689789239031448,
            0.81538923833917625, 0.833083879888400824, 0.849964527879591284, 0.866014688497164623,
            0.881218679385018416, 0.89556164497072699, 0.90902957098252969, 0.921609298145333953,
            0.933288535043079546, 0.94405587013625598, 0.953900782925491743, 0.96281365425581553,
            0.970785775763706332, 0.977809358486918289, 0.983877540706057016, 0.98898439524299175,
            0.99312493703744346, 0.99629513473312515, 0.998491950639595818, 0.99971372677344123};
}

std::vector<double> get_gl_grid(double cutoff){
    // GL grid shifted to (0, cutoff)
    auto grid = get_orig_gl_grid();
    double scale = cutoff/2;
    for(double & grid_val: grid){
        grid_val = scale * (grid_val + 1); // (b-a)/2 * x + (b+a)/2
    }
    return grid;
}
