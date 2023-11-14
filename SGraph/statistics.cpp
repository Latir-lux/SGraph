// 运行指令为./statistics [path] 例：./statistics ./pathlog/log/
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include <filesystem>

const int BLOCK_SIZE = 20000; // 定义块的大小

void processFile(const std::string& fileName, std::unordered_map<int, int>& blockCounts, int& totalVisits) {
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open input file " << fileName << std::endl;
        return;
    }

    int number;
    while (inputFile >> number) {
        int block = (number - 1) / BLOCK_SIZE; // 计算所属的块
        blockCounts[block]++;
        totalVisits++;
    }

    inputFile.close();
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_folder>" << std::endl;
        return 1;
    }

    std::string inputFolder = argv[1];

    std::unordered_map<int, int> blockCounts; // 用于统计每个块的访问次数
    int totalVisits = 0; // 记录总的访问次数

    for (const auto& entry : std::filesystem::directory_iterator(inputFolder)) {
        if (entry.path().extension() == ".txt") {
            processFile(entry.path(), blockCounts, totalVisits);
        }
    }

    // 创建输出文件
    std::string outputFileName = inputFolder + "/statistics.txt";
    std::ofstream outputFile(outputFileName);

    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to create output file." << std::endl;
        return 1;
    }

    double totalProbability = 0.0;
    double maxProbability = 0.0;

    for (const auto& pair : blockCounts) {
        double probability = static_cast<double>(pair.second) / totalVisits * 100.0; // 将概率转为百分比形式
        totalProbability += probability;
        if (probability > maxProbability) {
            maxProbability = probability;
        }

        // 输出时将概率乘以100
        outputFile << "Block " << pair.first << ": " << std::fixed << std::setprecision(2) << probability << "%" << std::endl;
    }

    double averageProbability = totalProbability / blockCounts.size();

    // 输出时将概率乘以100
    outputFile << "Average Probability: " << std::fixed << std::setprecision(2) << averageProbability << "%" << std::endl;
    outputFile << "Max Probability: " << std::fixed << std::setprecision(2) << maxProbability << "%" << std::endl;

    outputFile.close();

    std::cout << "Statistics saved to " << outputFileName << std::endl;

    return 0;
}
