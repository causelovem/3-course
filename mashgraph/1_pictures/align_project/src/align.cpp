#include "align.h"
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <tgmath.h>

#define PI 3.1415926535897932384626433832795
#define EPS 0.001 

using std::vector;
using std::tuple;
using std::get;
using std::tie;
using std::make_tuple;
using std::string;
using std::cout;
using std::endl;

class UnaryMatrixOpTuple /*+++++*/
{
    public:
        Matrix<double> kernel;
        uint radius;

        UnaryMatrixOpTuple(Matrix<double> kern):
        kernel{kern}, radius{(kern.n_rows - 1) / 2}
        {};

        tuple<uint, uint, uint> operator () (const Image &img)
        {
            uint size = 2 * radius + 1;
            uint r = 0, g = 0, b = 0;
            double sr = 0, sg = 0, sb = 0;

            for (uint i = 0; i < size; i++)
                for (uint j = 0; j < size; j++)
                {
                    tie(r, g, b) = img(i, j);
                    sr += r * kernel(i, j);
                    sg += g * kernel(i, j);
                    sb += b * kernel(i, j);
                }

            if (sr < 0)
                sr = 0;
            if (sr > 255)
                sr = 255;

            if (sg < 0)
                sg = 0;
            if (sg > 255)
                sg = 255;

            if (sb < 0)
                sb = 0;
            if (sb > 255)
                sb = 255;

            return make_tuple(static_cast<uint> (sr + 0.5), static_cast<uint> (sg + 0.5), static_cast<uint> (sb + 0.5));
        }
};

struct /*+++++*/
{
    /*http://www.cplusplus.com/reference/algorithm/sort/*/

    bool operator() (const std::tuple<uint, uint, uint> tup1, const std::tuple<uint, uint, uint> tup2)
    {
        double sum1 = (std::get<0>(tup1) + std::get<1>(tup1) + std::get<2>(tup1)) / 3;
        double sum2 = (std::get<0>(tup2) + std::get<1>(tup2) + std::get<2>(tup2)) / 3;

        if (sum1 < sum2)
            return true;
        else
            return false;
    }
} tuple_comp;

struct elem /*+++++*/
{
    double bright;
    uint num;

    bool operator == (const elem &el)
    {
        return fabs(el.bright - bright) < EPS;
    }
};

bool compare (const elem &el1, const elem &el2) /*+++++*/
{
    return el1.bright < el2.bright;
}

class UnaryMatrixOpMedian /*+++++*/
{
    public:
        int radius;
        UnaryMatrixOpMedian (int rad = 1) : radius{rad}
        {};

        tuple<uint, uint, uint> operator () (const Image &img)
        {
            int size = 2 * radius + 1;
            double r1 = 0, g1 = 0, b1 = 0;
            double r2 = 0, g2 = 0, b2 = 0;
            tuple<uint, uint, uint> tmp1, tmp2;
            std::vector<tuple<uint, uint, uint>> vec;

            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    vec.push_back(img(i, j));

            std::sort(vec.begin(), vec.end(), tuple_comp);

            tmp1 = vec[(vec.size() - 1) / 2];
            tmp2 = img(radius, radius);

            tie(r1, g1, b1) = tmp1;
            tie(r2, g2, b2) = tmp2;
            
            if (((r1 + g1 + b1) / 3) < ((r2 + g2 + b2) / 3))
                return tmp1;
            else
                return tmp2;
        }
};

Image mirroring(Image src_image, uint radius) /*+++++*/
{
    Image tmp_image(src_image.n_rows + 2 * radius, src_image.n_cols + 2 * radius);

    for (uint i = radius; i < src_image.n_rows + radius; i++)
        for (uint j = radius; j < src_image.n_cols + radius; j++)
            tmp_image(i, j) = src_image(i - radius, j - radius);

    for (uint i = 0; i < radius; i++)
        for (uint j = radius; j < tmp_image.n_cols - radius; j++)
        {
            tmp_image(i, j) = tmp_image(2 * radius - i, j);
            tmp_image(tmp_image.n_rows - 1 - i, j) = tmp_image(tmp_image.n_rows - 1 - 2 * radius + i, j);
        }
    
    for (uint j = 0; j < radius; j++)
        for (uint i = 0; i < tmp_image.n_rows; i++)
        {
            tmp_image(i, j) = tmp_image(i, 2 * radius - j);
            tmp_image(i, tmp_image.n_cols - 1 - j) = tmp_image(i, tmp_image.n_cols - 1 - 2 * radius + j);
        }

    return tmp_image;
}

void image_crop (Image &img_crop1, Image &img_crop2, int x_cord, int y_cord) /*+++++*/
{
    if ((x_cord >= 0) && (y_cord >= 0))
    {
        img_crop1 = img_crop1.submatrix(0, 0 + x_cord, img_crop1.n_rows - y_cord, img_crop1.n_cols - x_cord);
        img_crop2 = img_crop2.submatrix(0 + y_cord, 0, img_crop2.n_rows - y_cord, img_crop2.n_cols - x_cord);
    }
    else
    if ((x_cord >= 0) && (y_cord <= 0))
    {
        y_cord = abs(y_cord);
        img_crop1 = img_crop1.submatrix(0 + y_cord, 0 + x_cord, img_crop1.n_rows - y_cord, img_crop1.n_cols - x_cord);
        img_crop2 = img_crop2.submatrix(0, 0, img_crop2.n_rows - y_cord, img_crop2.n_cols - x_cord);
    }
    else
    if ((x_cord <= 0) && (y_cord <= 0))
    {
        x_cord = abs(x_cord);
        y_cord = abs(y_cord);
        img_crop1 = img_crop1.submatrix(0 + y_cord, 0, img_crop1.n_rows - y_cord, img_crop1.n_cols - x_cord);
        img_crop2 = img_crop2.submatrix(0, 0 + x_cord, img_crop2.n_rows - y_cord, img_crop2.n_cols - x_cord);
    }
    else
    if ((x_cord <= 0) && (y_cord >= 0))
    {
        x_cord = abs(x_cord);
        img_crop1 = img_crop1.submatrix(0, 0, img_crop1.n_rows - y_cord, img_crop1.n_cols - x_cord);
        img_crop2 = img_crop2.submatrix(0 + y_cord, 0 + x_cord, img_crop2.n_rows - y_cord, img_crop2.n_cols - x_cord);
    }
}

void find_metric1 (Image img_crop1, Image img_crop2, int &x_cord, int &y_cord) /*+++++*/
{
    uint min = -1, tmp_min = 0;
    Image tmp1, tmp2;

    for (int i = -15; i < 16; i++)
        for (int j = -15; j < 16; j++)
        {
            tmp1 = img_crop1;
            tmp2 = img_crop2;
            image_crop(tmp1, tmp2, i, j);

            for (uint k = 0; k < tmp1.n_rows; k++)
                for (uint l = 0; l < tmp1.n_cols; l++)
                    tmp_min += (((std::get<0>(tmp1(k,l)) + std::get<1>(tmp1(k,l)) + std::get<2>(tmp1(k,l))) - (std::get<0>(tmp2(k,l)) + std::get<1>(tmp2(k,l)) + std::get<2>(tmp2(k,l)))) 
                    * ((std::get<0>(tmp1(k,l)) + std::get<1>(tmp1(k,l)) + std::get<2>(tmp1(k,l))) - (std::get<0>(tmp2(k,l)) + std::get<1>(tmp2(k,l)) + std::get<2>(tmp2(k,l))))) / 9;

            tmp_min /= tmp1.n_cols * tmp1.n_rows;

            if (min > tmp_min)
            {
                min = tmp_min;
                x_cord = i;
                y_cord = j;
                //cout << min << endl;
            }

            tmp_min = 0;
        }

    //cout << min << endl;    

    return; 
}

void find_metric2 (Image img_crop1, Image img_crop2, int &x_cord, int &y_cord) /*+++++*/
{
    uint max = 0, tmp_max = 0;
    Image tmp1 = img_crop1, tmp2 = img_crop2;

    for (int i = -15; i < 16; i++)
        for (int j = -15; j < 16; j++)
        {
            tmp1 = img_crop1;
            tmp2 = img_crop2;
            image_crop(tmp1, tmp2, i, j);

            for (uint k = 0; k < tmp1.n_rows; k++)
                for (uint l = 0; l < tmp1.n_cols; l++)
                    tmp_max += (((std::get<0>(tmp1(k,l)) + std::get<1>(tmp1(k,l)) + std::get<2>(tmp1(k,l))) / 3) 
                    * ((std::get<0>(tmp2(k,l)) + std::get<1>(tmp2(k,l)) + std::get<2>(tmp2(k,l))) / 3));

            if (max < tmp_max)
            {
                max = tmp_max;
                x_cord = i;
                y_cord = j;
                //cout << max << endl;
            }

            tmp_max = 0;
        }

    //cout << max << endl;

    return; 
}

Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale) /*+++++*/
{
    int x_cord = 0, y_cord = 0;
    uint height = srcImage.n_rows / 3;
    Image img_red, img_blue, img_green, tmp;

    img_blue = srcImage.submatrix(0, 0, height, srcImage.n_cols);
    img_green = srcImage.submatrix(height, 0, height, srcImage.n_cols);
    img_red = srcImage.submatrix(2 * height, 0, height, srcImage.n_cols);

    double pers = 0.03;
    uint top_cut = img_blue.n_rows * pers, side_cut = img_blue.n_cols * pers;

    img_blue = img_blue.submatrix(top_cut, side_cut, img_blue.n_rows - 2 * top_cut, img_blue.n_cols - 2 * side_cut);
    img_green = img_green.submatrix(top_cut, side_cut, img_green.n_rows - 2 * top_cut, img_green.n_cols - 2 * side_cut);
    img_red = img_red.submatrix(top_cut, side_cut, img_red.n_rows - 2 * top_cut, img_red.n_cols - 2 * side_cut);

    tmp = img_blue;

    find_metric1(img_red, img_green, x_cord, y_cord);
    //find_metric2(img_red, img_green, x_cord, y_cord);

    image_crop(img_red, img_green, x_cord, y_cord);
    image_crop(img_blue, tmp, x_cord, y_cord);

    find_metric1(img_green, img_blue, x_cord, y_cord);
    //find_metric2(img_green, img_blue, x_cord, y_cord);

    for (uint k = 0; k < img_red.n_rows; k++)
        for (uint l = 0; l < img_red.n_cols; l++)
            std::get<0>(img_green(k,l)) = std::get<0>(img_red(k,l));

    image_crop(img_green, img_blue, x_cord, y_cord);

    for (uint k = 0; k < img_green.n_rows; k++)
        for (uint l = 0; l < img_green.n_cols; l++)
            std::get<2>(img_green(k,l)) = std::get<2>(img_blue(k,l));

    srcImage = img_green;
   
    if (isPostprocessing == true)
    {
        if (postprocessingType == "--gray-world")
            srcImage = gray_world(srcImage);

        if (postprocessingType == "--unsharp")
            srcImage = unsharp(srcImage);

        if (postprocessingType == "--autocontrast")
            srcImage = autocontrast(srcImage, fraction);
    }

    return srcImage;
}

Image sobel_x(Image src_image) /*+++++*/
{
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};

    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) /*+++++*/
{
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};

    return custom(src_image, kernel);
}

Image unsharp(Image src_image) /*+++++*/
{
    Matrix<double> kernel = {{-0.16666666667, -0.6666666667, -0.16666666667},
                             {-0.6666666667, 4.3333333333, -0.6666666667},
                             {-0.16666666667, -0.6666666667, -0.16666666667}};

    return custom(src_image, kernel);
}

Image gray_world(Image src_image) /*+++++*/
{
    double sr = 0, sg = 0, sb = 0, s = 0, tmp = 0;

    for (uint i = 0; i < src_image.n_rows; i++)
        for (uint j = 0; j < src_image.n_cols; j++)
        {
            sr += std::get<0>(src_image(i,j));
            sg += std::get<1>(src_image(i,j));
            sb += std::get<2>(src_image(i,j));
        }

    s = src_image.n_rows * src_image.n_cols;
    sr /= s;
    sg /= s;
    sb /= s;

    s = (sr + sg + sb) / 3;

    sr = s / sr;
    sg = s / sg;
    sb = s / sb;

    for (uint i = 0; i < src_image.n_rows; i++)
        for (uint j = 0; j < src_image.n_cols; j++)
        {
            if ((tmp = std::get<0>(src_image(i,j)) * sr) > 255)
                std::get<0>(src_image(i,j)) = 255;
            else
                std::get<0>(src_image(i,j)) = tmp;

            if ((tmp = std::get<1>(src_image(i,j)) * sg) > 255)
                std::get<1>(src_image(i,j)) = 255;
            else
                std::get<1>(src_image(i,j)) = tmp;

            if ((tmp = std::get<2>(src_image(i,j)) * sb) > 255)
                std::get<2>(src_image(i,j)) = 255;
            else
                std::get<2>(src_image(i,j)) = tmp;
        }

    return src_image;
}

Image resize(Image src_image, double scale)
{
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) /*+++++*/
{
    uint radius = (kernel.n_rows - 1) / 2;
    Image tmp_image(src_image.n_rows + 2 * radius, src_image.n_cols + 2 * radius);
    UnaryMatrixOpTuple op(kernel);

    tmp_image = mirroring(src_image, radius);
    tmp_image = tmp_image.unary_map(op);
    src_image = tmp_image.submatrix(radius, radius, src_image.n_rows, src_image.n_cols);

    return src_image;
}

Image autocontrast(Image src_image, double fraction) /*+++++*/
{
    std::vector<elem> vec;
    std::vector<elem>::iterator pos;
    elem tmp_el;
    double tmp1 = 0, tmp2 = 0, min = 0, max = 0;
    uint r = 0, g = 0, b = 0;
    uint cut = src_image.n_rows * src_image.n_cols * fraction;

    for (uint i = 0; i < src_image.n_rows; i++)
        for (uint j = 0; j < src_image.n_cols; j++)
        {
            tmp_el.num = 1;
            tie(r, g, b) = src_image(i, j);
            tmp_el.bright = 0.2125 * r + 0.7154 * g + 0.0721 * b;
            if ((pos = find(vec.begin(), vec.end(), tmp_el)) != vec.end())
                pos->num++;
            else
               vec.push_back(tmp_el); 
        }
    std::sort(vec.begin(), vec.end(), compare);

    for (uint i = 0; i < vec.size(); i++)
    {
        tmp1 += vec[i].num;
        tmp2 += vec[vec.size() - 1 - i].num;

        if ((min < EPS) && (tmp1 > cut))
            min = vec[i].bright;

        if ((max < EPS) && (tmp2 > cut))
            max = vec[vec.size() - 1 - i].bright;
    }   

    double consta = 255 / (max - min), nr = 0, ng = 0, nb = 0;
    for (uint i = 0; i < src_image.n_rows; i++)
        for (uint j = 0; j < src_image.n_cols; j++)
        {
            tie(r, g, b) = src_image(i, j);
            tmp1 = 0.2125 * r + 0.7154 * g + 0.0721 * b;

            nr = r * (1 - min / tmp1) * consta;
            ng = g * (1 - min / tmp1) * consta;
            nb = b * (1 - min / tmp1) * consta;

            if (nr < 0)
                nr = 0;
            if (nr > 255)
                nr = 255;

            if (ng < 0)
                ng = 0;
            if (ng > 255)
                ng = 255;

            if (nb < 0)
                nb = 0;
            if (nb > 255)
                nb = 255;

            src_image(i, j) = make_tuple(static_cast<uint> (nr + 0.5), static_cast<uint> (ng + 0.5), static_cast<uint> (nb + 0.5));
        }

    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius) /*+++++*/
{
    Matrix<double> gauss_kernel(2 * radius + 1, 2 * radius + 1);
    long double coeff = 2 * PI * sigma * sigma, tmp;
    long double norm = 0;

    for (uint i = 0; i < gauss_kernel.n_rows; i++)
        for (uint j = 0; j < gauss_kernel.n_cols; j++)
        {
            tmp = ((i - radius) * (i - radius) + (j - radius) * (j - radius)) / (2 * sigma * sigma);
            tmp = std::exp(-tmp);
            gauss_kernel(i, j) = tmp / coeff;
            norm += gauss_kernel(i, j);
        }

    for (uint i = 0; i < gauss_kernel.n_rows; i++)
        for (uint j = 0; j < gauss_kernel.n_cols; j++)
            gauss_kernel(i, j) /= norm;

    return custom(src_image, gauss_kernel);
}

Image gaussian_separable(Image src_image, double sigma, int radius) /*+++++*/
{
    Matrix<double> gauss_kernel(2 * radius + 1, 2 * radius + 1);
    Matrix<double> gauss_kernel_col(2 * radius + 1, 1), gauss_kernel_row(1, 2 * radius + 1);
    long double coeff = 2 * PI * sigma * sigma, tmp = 0;
    long double norm = 0;

    for (uint i = 0; i < gauss_kernel.n_rows; i++)
        for (uint j = 0; j < gauss_kernel.n_cols; j++)
        {
            tmp = ((i - radius) * (i - radius) + (j - radius) * (j - radius)) / (2 * sigma * sigma);
            tmp = std::exp(-tmp);
            gauss_kernel(i, j) = tmp / coeff;
            norm += gauss_kernel(i, j);
        }

    for (uint i = 0; i < gauss_kernel.n_rows; i++)
        for (uint j = 0; j < gauss_kernel.n_cols; j++)
            gauss_kernel(i, j) /= norm;

    for (uint i = 0; i < gauss_kernel.n_rows; i++)
    {
        gauss_kernel_col(i, 0) = 0;
        for (uint j = 0; j < gauss_kernel.n_cols; j++)
            gauss_kernel_col(i, 0) += gauss_kernel(i, j);
    }

    for (uint i = 0; i < gauss_kernel_col.n_rows; i++)
        gauss_kernel_row(0, i) = gauss_kernel_col(i, 0);

    Image tmp_image(src_image.n_rows + 2 * radius, src_image.n_cols + 2 * radius);
    Image res_image1(src_image.n_rows, src_image.n_cols);

    tmp_image = mirroring(src_image, radius);
    
    uint r = 0, g = 0, b = 0;
    double sr = 0, sg = 0, sb = 0;

    for (uint i = radius; i < tmp_image.n_rows - radius; i++)
        for (uint j = radius; j < tmp_image.n_cols - radius; j++)
        {
            r = g = b = 0;
            sr = sg = sb = 0;

            for (int k = -radius; k < radius + 1; k++)
            {
                tie(r, g, b) = tmp_image(i + k, j);
                sr += r * gauss_kernel_col(radius + k, 0);
                sg += g * gauss_kernel_col(radius + k, 0);
                sb += b * gauss_kernel_col(radius + k, 0);
            }

            if (sr < 0)
                sr = 0;
            if (sr > 255)
                sr = 255;

            if (sg < 0)
                sg = 0;
            if (sg > 255)
                sg = 255;

            if (sb < 0)
                sb = 0;
            if (sb > 255)
                sb = 255;

            res_image1(i - radius, j -radius) = make_tuple(static_cast<uint> (sr + 0.5), static_cast<uint> (sg + 0.5), static_cast<uint> (sb + 0.5));
        }

    tmp_image = mirroring(res_image1, radius);
    Image res_image(src_image.n_rows, src_image.n_cols);

    for (uint j = radius; j < tmp_image.n_cols - radius; j++)
        for (uint i = radius; i < tmp_image.n_rows - radius; i++)
        {
            r = g = b = 0;
            sr = sg = sb = 0;

            for (int k = -radius; k < radius + 1; k++)
            {
                tie(r, g, b) = tmp_image(i, j + k);
                sr += r * gauss_kernel_row(0, radius + k);
                sg += g * gauss_kernel_row(0, radius + k);
                sb += b * gauss_kernel_row(0, radius + k);
            }

            if (sr < 0)
                sr = 0;
            if (sr > 255)
                sr = 255;

            if (sg < 0)
                sg = 0;
            if (sg > 255)
                sg = 255;

            if (sb < 0)
                sb = 0;
            if (sb > 255)
                sb = 255;

            res_image(i - radius, j - radius) = make_tuple(static_cast<uint> (sr + 0.5), static_cast<uint> (sg + 0.5), static_cast<uint> (sb + 0.5));
        }

    return res_image;
}

Image median(Image src_image, int radius) /*+++++*/
{
    Image tmp_image(src_image.n_rows + 2 * radius, src_image.n_cols + 2 * radius);
    UnaryMatrixOpMedian op(radius);

    tmp_image = mirroring(src_image, radius);
    tmp_image = tmp_image.unary_map(op);
    src_image = tmp_image.submatrix(radius, radius, src_image.n_rows, src_image.n_cols);

    return src_image;
}

Image median_linear(Image src_image, int radius)
{
    return src_image;
}

Image median_const(Image src_image, int radius)
{
    return src_image;
}

Image canny(Image src_image, int threshold1, int threshold2) /*+++++*/
{
    /*Algorithm was partially taken from
    https://habrahabr.ru/post/114589/*/

    Image src_image_save = src_image.deep_copy();

    /*GAUSS*/
    //src_image = gaussian(src_image, 1.4, 2);
    src_image = gaussian_separable(src_image, 1.4, 2);

    /*SOBEL X Y*/
    Image Ix = sobel_x(src_image);
    Image Iy = sobel_y(src_image);
    Matrix<double> modG(Ix.n_rows, Ix.n_cols);
    Matrix<double> theta(Ix.n_rows, Ix.n_cols);

    for (uint i = 0; i < Ix.n_rows; i++)
        for (uint j = 0; j < Ix.n_cols; j++)
        {
            modG(i, j) = sqrt (std::get<0>(Ix(i, j)) * std::get<0>(Ix(i, j)) + std::get<0>(Iy(i, j)) * std::get<0>(Iy(i, j)));
            theta(i, j) = (-PI / 2) + (PI / 4) * round((4 / PI) * atan2(std::get<0>(Iy(i, j)), std::get<0>(Ix(i, j))));
        }

    /*NOT MAX*/
    for (uint i = 1; i < Ix.n_rows - 1; i++)
        for (uint j = 1; j < Ix.n_cols - 1; j++)
        {
            if (((abs(theta(i, j)) < EPS) || (abs(theta(i, j) - PI) < EPS) || (abs(theta(i, j) + PI) < EPS))
                && ((modG(i, j) < modG(i, j + 1)) || (modG(i, j) < modG(i, j - 1))))
                modG(i, j) = 0.0;
            else
            if (((abs(theta(i, j) - PI / 4) < EPS) || (abs(theta(i, j) + (3 / 4) * PI) < EPS))
                && ((modG(i, j) < modG(i - 1, j + 1)) || (modG(i, j) < modG(i + 1, j - 1))))
                modG(i, j) = 0.0;
            else
            if (((abs(theta(i, j) - PI / 2) < EPS) || (abs(theta(i, j) + PI / 2) < EPS))
                && (modG(i, j) < modG(i + 1, j) || (modG(i, j) < modG(i - 1, j))))
                modG(i, j) = 0.0;
            else
            if (((abs(theta(i, j) - (3 / 4) * PI) < EPS) || (abs(theta(i, j) + PI / 4) < EPS))
                && (modG(i, j) < modG(i - 1, j - 1) || (modG(i, j) < modG(i + 1, j + 1))))
                modG(i, j) = 0.0;
        }

    /*TWO THRESHOLD*/
    Matrix<uint> border_map(Ix.n_rows, Ix.n_cols);

    for (uint i = 1; i < Ix.n_rows - 1; i++)
        for (uint j = 1; j < Ix.n_cols - 1; j++)
        {
            if (modG(i, j) > threshold2)
                border_map(i, j) = 255;
            else
            if (modG(i, j) < threshold1)
                border_map(i, j) = 0;
            else
                border_map(i, j) = 127;
        }

    /*FIND BORDERS*/
    std::vector<tuple<uint, uint>> white_points;
    std::vector<tuple<uint, uint>> gray_points;
    Matrix<uint> border_map_new(Ix.n_rows, Ix.n_cols);

    for (uint i = 0; i < Ix.n_rows; i++)
        for (uint j = 0; j < Ix.n_cols; j++)
        {
            if (border_map(i, j) == 255)
                white_points.push_back(make_tuple(i, j));
            if (border_map(i, j) == 127)
                gray_points.push_back(make_tuple(i, j));
        }

    int i1 = 0, i2 = 0, j1 = 0, j2 = 0;
    for (uint i = 0; i < white_points.size(); i++)
        for (uint j = 0; j < gray_points.size(); j++)
        {
            tie(i1, j1) = white_points[i];
            tie(i2, j2) = gray_points[j];

            if ((abs(i1 - i2) <= 1) && (abs(j1 - j2) <= 1))
            {
                white_points.push_back(make_tuple(i2, j2));
                gray_points.erase(gray_points.begin() + j);
            }
        }

    for (uint i = 0; i < white_points.size(); i++)
    {
        tie(i1, j1) = white_points[i];
        border_map_new(i1, j1) = 255;
    }


    /*OUT BORDERS MAP*/
    for (uint i = 0; i < border_map_new.n_rows; i++)
        for (uint j = 0; j < border_map_new.n_cols; j++)
            src_image(i, j) = make_tuple(border_map_new(i, j),border_map_new(i, j),border_map_new(i, j));
    return src_image;

    /*FIND BORDER*/
    double pers = 0.05;
    uint top_border = border_map_new.n_rows * pers, side_border = border_map_new.n_cols * pers;
    uint top = 0, bottom = 0, left = 0, right = 0;
    uint max1 = 0, max2 = 0, max1_num = 0, max2_num = 0, tmp = 0;

    for (uint i = 0; i < top_border; i++)
    {
        for (uint j = 0; j < border_map_new.n_cols; j++)
            if (border_map_new(i, j) == 255)
                tmp++;

        if (tmp >= max2)
        {
            if (tmp >= max1)
            {
                max2 = max1;
                max2_num = max1_num;
                max1 = tmp;
                max1_num = i;
            }
            else
            {
                max2 = tmp;
                max2_num = i;
            }

            i++;
        }

        tmp = 0;
    }
    top = fmax(max1_num, max2_num);

    max1 = max2 = tmp = 0;
    max1_num = max2_num = border_map_new.n_rows - 1;

    for (uint i = border_map_new.n_rows - 1; i > border_map_new.n_rows - 1 - top_border; i--)
    {
        for (uint j = 0; j < border_map_new.n_cols; j++)
            if (border_map_new(i, j) == 255)
                tmp++;

        if (tmp >= max2)
        {
            if (tmp >= max1)
            {
                max2 = max1;
                max2_num = max1_num;
                max1 = tmp;
                max1_num = i;
            }
            else
            {
                max2 = tmp;
                max2_num = i;
            }

            i--;
        }

        tmp = 0;
    }
    bottom = fmin(max1_num, max2_num);

    max1 = max2 = max1_num = max2_num = tmp = 0;

    for (uint i = 0; i < side_border; i++)
    {
        for (uint j = 0; j < border_map_new.n_rows; j++)
            if (border_map_new(j, i) == 255)
                tmp++;

        if (tmp >= max2)
        {
            if (tmp >= max1)
            {
                max2 = max1;
                max2_num = max1_num;
                max1 = tmp;
                max1_num = i;
            }
            else
            {
                max2 = tmp;
                max2_num = i;
            }

            i++;
        }

        tmp = 0;
    }
    left = fmax(max1_num, max2_num);

    max1 = max2 = tmp = 0;
    max1_num = max2_num = border_map_new.n_cols - 1;

    for (uint i = border_map_new.n_cols - 1; i > border_map_new.n_cols - 1 - side_border; i--)
    {
        for (uint j = 0; j < border_map_new.n_rows; j++)
            if (border_map_new(j, i) == 255)
                tmp++;

        if (tmp >= max2)
        {
            if (tmp >= max1)
            {
                max2 = max1;
                max2_num = max1_num;
                max1 = tmp;
                max1_num = i;
            }
            else
            {
                max2 = tmp;
                max2_num = i;
            }

            i--;
        }

        tmp = 0;
    }
    right = fmin(max1_num, max2_num);

    /*CUT*/
    src_image_save = src_image_save.submatrix(top, left, bottom - top, right - left);

    return src_image_save;
}