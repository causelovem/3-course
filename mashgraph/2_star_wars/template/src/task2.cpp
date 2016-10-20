#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cmath>

#include "classifier.h"
#include "EasyBMP.h"
#include "linear.h"
#include "argvparser.h"
#include "matrix.h"

#define PI 3.1415926535897932384626433832795
#define EPS 0.001
#define SEG_NUM 16
#define BLOCK_SIZE 10

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;

using CommandLineProcessing::ArgvParser;

typedef vector<pair<BMP*, int> > TDataSet;
typedef vector<pair<string, int> > TFileList;
typedef vector<pair<vector<float>, int> > TFeatures;

class UnaryMatrixOp /*+++++*/
{
    public:
        Matrix <float> kernel;
        size_t vert_radius;
        size_t hor_radius;

        UnaryMatrixOp(Matrix<float> kern):
        kernel{kern}, vert_radius{(kern.n_rows - 1) / 2}, hor_radius{(kern.n_cols - 1) / 2}
        {};

        float operator () (const Matrix <float> &img)
        {
            float sum = 0.0;

            for (size_t i = 0; i < img.n_rows; i++)
                for (size_t j = 0; j < img.n_cols; j++) 
                    sum += img(i, j) * kernel(i, j);

            return sum;
        }
};

// Load list of files and its labels from 'data_file' and
// stores it in 'file_list'
void LoadFileList(const string& data_file, TFileList* file_list)
{
    ifstream stream(data_file.c_str());

    string filename;
    int label;
    
    int char_idx = data_file.size() - 1;
    for (; char_idx >= 0; --char_idx)
        if (data_file[char_idx] == '/' || data_file[char_idx] == '\\')
            break;
    string data_path = data_file.substr(0,char_idx+1);
    
    while(!stream.eof() && !stream.fail()) {
        stream >> filename >> label;
        if (filename.size())
            file_list->push_back(make_pair(data_path + filename, label));
    }

    stream.close();
}

// Load images by list of files 'file_list' and store them in 'data_set'
void LoadImages(const TFileList& file_list, TDataSet* data_set)
{
    for (size_t img_idx = 0; img_idx < file_list.size(); ++img_idx) {
            // Create image
        BMP* image = new BMP();
            // Read image from file
        image->ReadFromFile(file_list[img_idx].first.c_str());
            // Add image and it's label to dataset
        data_set->push_back(make_pair(image, file_list[img_idx].second));
    }
}

// Save result of prediction to file
void SavePredictions(const TFileList& file_list,
                     const TLabels& labels, 
                     const string& prediction_file)
{
        // Check that list of files and list of labels has equal size 
    assert(file_list.size() == labels.size());
        // Open 'prediction_file' for writing
    ofstream stream(prediction_file.c_str());

        // Write file names and labels to stream
    for (size_t image_idx = 0; image_idx < file_list.size(); ++image_idx)
        stream << file_list[image_idx].first << " " << labels[image_idx] << endl;
    stream.close();
}

void LBP (std::vector<float> &descriptor, const Matrix <float> &bright, int block_rows, int block_cols) /*+++++*/
{
    Matrix <float> new_bright = bright.deep_copy();

    new_bright = new_bright.extra_borders(1, 1);

    /*CELL COMPUTUNG*/
    std::vector<std::vector<float>> tmp_descriptor;
    for (uint i = 1; i < new_bright.n_rows - 1; i += block_rows)
    {
        for (uint j = 1; j < new_bright.n_cols - 1; j += block_cols)
        {
            std::vector<float> one_image_features;
            one_image_features.resize(256, 0);

            for (uint k = i; k < i + block_rows; k++)
            {
                for (uint l = j; l < j + block_cols; l++)
                {
                    int num = 0;

                    num += (1 << 7) * (new_bright(k - 1, l - 1) >= new_bright(k, l));
                    num += (1 << 6) * (new_bright(k - 1, l - 0) >= new_bright(k, l));
                    num += (1 << 5) * (new_bright(k - 1, l + 1) >= new_bright(k, l));
                    num += (1 << 4) * (new_bright(k - 0, l + 1) >= new_bright(k, l));
                    num += (1 << 3) * (new_bright(k + 1, l + 1) >= new_bright(k, l));
                    num += (1 << 2) * (new_bright(k + 1, l - 0) >= new_bright(k, l));
                    num += (1 << 1) * (new_bright(k + 1, l - 1) >= new_bright(k, l));
                    num += (1 << 0) * (new_bright(k - 0, l - 1) >= new_bright(k, l));
                    one_image_features[num]++;
                }
            }
            tmp_descriptor.push_back(one_image_features);
        }
    }

    /*NORM*/
    size_t size = tmp_descriptor.size();
    for (size_t i = 0; i < size; i++)
    {
        float sum = 0.0;
        for (size_t j = 0; j < tmp_descriptor[i].size(); j++)
            sum += tmp_descriptor[i][j] * tmp_descriptor[i][j];
        
        if ((sum = sqrt(sum)) > EPS)
            for (size_t j = 0; j < tmp_descriptor[i].size(); j++)
                tmp_descriptor[i][j] /= sum;
    }

    /*MERGE*/
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < tmp_descriptor[i].size(); j++)
            descriptor.push_back(tmp_descriptor[i][j]);

    return;
}

void color_feature (std::vector<float> &descriptor, BMP *img) /*+++++*/
{
    int height = img->TellHeight(), width = img->TellWidth();
    int block_rows = height / 8, block_cols = width / 8;
    int cut_rows = height % 8, cut_cols = width % 8;

    /*CELL COMPUTUNG*/
    std::vector<std::vector<float>> tmp_descriptor;
    for (int i = cut_rows / 2; i < height - cut_rows; i += block_rows)
    {
        for (int j = cut_cols / 2; j < width - cut_cols; j += block_cols)
        {
            std::vector<float> one_image_features;
            one_image_features.resize(3, 0);

            for (int k = i; k < i + block_rows; k++)
            {
                for (int l = j; l < j + block_cols; l++)
                {
                    one_image_features[0] += img->GetPixel(l, k).Red;
                    one_image_features[1] += img->GetPixel(l, k).Green;
                    one_image_features[2] += img->GetPixel(l, k).Blue;
                }
            }
            one_image_features[0] /= 64;
            one_image_features[1] /= 64;
            one_image_features[2] /= 64;

            tmp_descriptor.push_back(one_image_features);
        }
    }

    /*NORM*/
    size_t size = tmp_descriptor.size();
    for (size_t i = 0; i < size; i++)      
        for (size_t j = 0; j < tmp_descriptor[i].size(); j++)
            tmp_descriptor[i][j] /= 255;

    /*MERGE*/
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < tmp_descriptor[i].size(); j++)
            descriptor.push_back(tmp_descriptor[i][j]);

    return;
}

// Exatract features from dataset.
// You should implement this function by yourself =)
void ExtractFeatures(const TDataSet& data_set, TFeatures* features) /*+++++*/
{
    Matrix <float> sobel_x(1,3);
    Matrix <float> sobel_y(3,1);

    sobel_x(0, 0) = -1.0;
    sobel_x(0, 1) = 0;
    sobel_x(0, 2) = 1.0;

    sobel_y(0, 0) = -1.0;
    sobel_y(1, 0) = 0;
    sobel_y(2, 0) = 1.0;

    for (size_t image_idx = 0; image_idx < data_set.size(); image_idx++)
    {
        BMP *img;
        img = std::get<0>(data_set[image_idx]);
        int height = img->TellHeight(), width = img->TellWidth();
        Matrix <float> bright(height, width);

        for (int i = 0; i < width; i++)
            for (int j = 0; j < height; j++)
                bright(j, i) = 0.299 * img->GetPixel(i, j).Red + 0.587 * img->GetPixel(i, j).Green + 0.114 * img->GetPixel(i, j).Blue;

        Matrix <float> Ix(height, width);
        Matrix <float> Iy(height, width);

        UnaryMatrixOp opx(sobel_x);
        UnaryMatrixOp opy(sobel_y);
        Ix = bright.unary_map(opx);
        Iy = bright.unary_map(opy);

        Matrix <float> modG(height, width);
        Matrix <float> theta(height, width);

        for (size_t i = 0; i < Ix.n_rows; i++)
            for (size_t j = 0; j < Ix.n_cols; j++)
            {
                modG(i, j) = sqrt (Ix(i, j) * Ix(i, j) + Iy(i, j) * Iy(i, j));
                theta(i, j) = atan2f(Iy(i, j), Ix(i, j));
            }

        int block_rows = modG.n_rows / BLOCK_SIZE, block_cols = modG.n_cols / BLOCK_SIZE;
        int cut_rows = modG.n_rows % BLOCK_SIZE, cut_cols = modG.n_cols % BLOCK_SIZE;

        modG = modG.submatrix(cut_rows / 2, cut_cols / 2, modG.n_rows - cut_rows, modG.n_cols - cut_cols);
        theta = theta.submatrix(cut_rows / 2, cut_cols / 2, theta.n_rows - cut_rows, theta.n_cols - cut_cols);
        bright = bright.submatrix(cut_rows / 2, cut_cols / 2, bright.n_rows - cut_rows, bright.n_cols - cut_cols);

        /*CELL COMPUTUNG*/
        std::vector<std::vector<float>> tmp_descriptor;
        for (uint i = 0; i < theta.n_rows; i += block_rows)
        {
            for (uint j = 0; j < theta.n_cols; j += block_cols)
            {
                std::vector<float> one_image_features;
                one_image_features.resize(SEG_NUM, 0);

                for (uint k = i; k < i + block_rows; k++)
                {
                    for (uint l = j; l < j + block_cols; l++)
                    {
                        int seg = fabsf(theta(k, l)) * (SEG_NUM / 2) / PI;
                        if (theta(k, l) < 0)
                        {
                            if (fabs(fabsf(theta(k, l)) - PI) < EPS)
                                one_image_features[0] += modG(k, l);
                            else
                                one_image_features[(SEG_NUM / 2) - 1 - seg] += modG(k, l);
                        }
                        else
                        if (theta(k, l) >= 0)
                        {
                            if (fabs(fabsf(theta(k, l)) - PI) < EPS)
                                one_image_features[one_image_features.size() - 1] += modG(k, l);
                            else
                                one_image_features[(SEG_NUM / 2) + seg] += modG(k, l);
                        }
                    }
                }
                tmp_descriptor.push_back(one_image_features);
            }
        }

        /*NORM*/
        size_t size = tmp_descriptor.size();
        for (size_t i = 0; i < size; i++)
        {
            float sum = 0.0;
            for (size_t j = 0; j < tmp_descriptor[i].size(); j++)
                sum += tmp_descriptor[i][j] * tmp_descriptor[i][j];
            
            if ((sum = sqrt(sum)) > EPS)
                for (size_t j = 0; j < tmp_descriptor[i].size(); j++)
                    tmp_descriptor[i][j] /= sum;
        }

        /*MERGE*/
        std::vector<float> descriptor;
        for (size_t i = 0; i < size; i++)
            for (size_t j = 0; j < tmp_descriptor[i].size(); j++)
                descriptor.push_back(tmp_descriptor[i][j]);

        LBP(descriptor, bright, block_rows, block_cols);
        color_feature(descriptor, img);

        features->push_back(make_pair(descriptor, std::get<1>(data_set[image_idx])));
    }

    return;
}

// Clear dataset structure
void ClearDataset(TDataSet* data_set)
{
        // Delete all images from dataset
    for (size_t image_idx = 0; image_idx < data_set->size(); ++image_idx)
        delete (*data_set)[image_idx].first;
        // Clear dataset
    data_set->clear();
}

// Train SVM classifier using data from 'data_file' and save trained model
// to 'model_file'
void TrainClassifier(const string& data_file, const string& model_file)
{
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // Model which would be trained
    TModel model;
        // Parameters of classifier
    TClassifierParams params;
    
        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // PLACE YOUR CODE HERE
        // You can change parameters of classifier here
    params.C = 0.01;
    TClassifier classifier(params);
        // Train classifier
    classifier.Train(features, &model);
        // Save model to file
    model.Save(model_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

// Predict data from 'data_file' using model from 'model_file' and
// save predictions to 'prediction_file'
void PredictData(const string& data_file,
                 const string& model_file,
                 const string& prediction_file)
{
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // List of image labels
    TLabels labels;

        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // Classifier 
    TClassifier classifier = TClassifier(TClassifierParams());
        // Trained model
    TModel model;
        // Load model from file
    model.Load(model_file);
        // Predict images by its features using 'model' and store predictions
        // to 'labels'
    classifier.Predict(features, model, &labels);

        // Save predictions
    SavePredictions(file_list, labels, prediction_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

int main(int argc, char** argv)
{
    // Command line options parser
    ArgvParser cmd;
        // Description of program
    cmd.setIntroductoryDescription("Machine graphics course, task 2. CMC MSU, 2014.");
        // Add help option
    cmd.setHelpOption("h", "help", "Print this help message");
        // Add other options
    cmd.defineOption("data_set", "File with dataset",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("model", "Path to file to save or load model",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("predicted_labels", "Path to file to save prediction results",
        ArgvParser::OptionRequiresValue);
    cmd.defineOption("train", "Train classifier");
    cmd.defineOption("predict", "Predict dataset");
        
        // Add options aliases
    cmd.defineOptionAlternative("data_set", "d");
    cmd.defineOptionAlternative("model", "m");
    cmd.defineOptionAlternative("predicted_labels", "l");
    cmd.defineOptionAlternative("train", "t");
    cmd.defineOptionAlternative("predict", "p");

        // Parse options
    int result = cmd.parse(argc, argv);

        // Check for errors or help option
    if (result) {
        cout << cmd.parseErrorDescription(result) << endl;
        return result;
    }

        // Get values 
    string data_file = cmd.optionValue("data_set");
    string model_file = cmd.optionValue("model");
    bool train = cmd.foundOption("train");
    bool predict = cmd.foundOption("predict");

        // If we need to train classifier
    if (train)
        TrainClassifier(data_file, model_file);
        // If we need to predict data
    if (predict) {
            // You must declare file to save images
        if (!cmd.foundOption("predicted_labels")) {
            cerr << "Error! Option --predicted_labels not found!" << endl;
            return 1;
        }
            // File to save predictions
        string prediction_file = cmd.optionValue("predicted_labels");
            // Predict data
        PredictData(data_file, model_file, prediction_file);
    }
}