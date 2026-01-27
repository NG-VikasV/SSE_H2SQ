#ifndef RESULT_WRITER_HPP
#define RESULT_WRITER_HPP

#include <string>
#include <vector>

class ResultWriter {
public:
    // Constructor: choose output file
    explicit ResultWriter(const std::string& filename,
                          bool append = false);

    // Write a single row of results
    void write_row(const std::vector<std::string>& headers,
                   const std::vector<double>& values);

    // Optional: write comment or metadata line
    void write_comment(const std::string& comment);

private:
    std::string filename_;
    bool append_;
};

#endif // RESULT_WRITER_HPP
