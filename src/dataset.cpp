#include "dataset.h"

Dataset::Dataset(const std::string &dataset_name) : name_(dataset_name)
{
    static const unordered_map<std::string, std::string> name_to_filename =
        {
            // {"Adult", "../datasets/adultSky.txt"},
            {"Adult", "../datasets/Adult.txt"},
            {"Anti-Cor_4_10000", "../datasets/Anti-Cor_4_10000.txt"},
            {"German", "../datasets/German.txt"},
            {"Compas", "../datasets/Compas.txt"},
            {"smart-home", "../datasets/smart-home.txt"},
            {"smart-home-skyline", "../datasets/smart-home-skyline.txt"},
            {"healthcare", "../datasets/healthcare.txt"},
            {"wearable", "../datasets/wearable.txt"}
            };


    if (!name_to_filename.count(name_))
    {
        Fail("Unknown dataset name: " + name_);
    }
    const std::string &file_name = name_to_filename.at(name_);

    std::cerr << "Reading dataset from " << file_name << "..." << std::endl;
    std::ifstream fin;
    fin.open(file_name,std::ios::in);
    if (!fin)
    {
        std::cout << "Cannot open file " << file_name << " for reading \n";
        exit(1);
    }
    std::string first_line;
    getline(fin,first_line);
    std::stringstream first_line_stream(first_line);
    size_t D_size;
    first_line_stream >> D_size;
    D_size_ = D_size;
    std::cout << "D_size:" << D_size << std::endl;
    this->num_points_ = 1;
    size_t dim;
    first_line_stream >> dim;
    dim_ = dim;
    std::cout << "dim:" << dim << std::endl;
    // size_t attribute_num;
    // first_line_stream >> attribute_num;
    // att_num_ = attribute_num;
    // std::cout << "att_num:" << attribute_num << std::endl;

    std::string line;
    while (std::getline(fin, line))
    {
        if (line.empty())
            continue;
        std::istringstream iss(line);
        vector<double> coordinates(dim);
        // vector<std::string> attributes(attribute_num);

        for (size_t i = 0; i < dim; ++i)
        {
            if (!(iss >> coordinates[i]))
            {
                std::cerr << "Error reading point coordinates" << std::endl;
                exit(1);
            }
        }
        // for (size_t i = 0; i < attribute_num; ++i)
        // {
        //     if (!(iss >> attributes[i]))
        //     {
        //         std::cerr << "Error reading point attributes" << std::endl;
        //         exit(1);
        //     }
        // }
        // Point p((num_points_-1), attribute_num, dim, attributes, coordinates);
        Point p((num_points_-1), dim,  coordinates);
        points_.emplace_back(p);

        num_points_++;
    }

    // int j = 0;
    // for (const auto& point : points_) {
    //     std::cout << "Point: " ;
    //     for (size_t i = 0; i < dim_; ++i) {
    //         std::cout << points_[j].get_coordinate(i);
    //         if (i < dim_ - 1) {
    //             std::cout << ", "; // Add a comma between coordinates
    //         }
    //     }
    //     j = j+1;
    //     std::cout<< std::endl;
    // }

    std::cerr << "Already read dataset " << name_ << " with " << num_points_-- << " points." << std::endl;
}

const vector<Point> &Dataset::get_universe_points() const
{
    return points_;
}

const std::string &Dataset::get_name() const
{
    return name_;
}
