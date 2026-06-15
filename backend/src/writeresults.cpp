#include "../include/writeresults.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// --------------------------------------------------
// Constructor
// --------------------------------------------------
ResultWriter::ResultWriter(const std::string& filename, bool append)
    : filename_(filename), append_(append)
{}

// --------------------------------------------------
// Write one line of observables (append-safe)
// --------------------------------------------------
void ResultWriter::write_row(const std::vector<std::string>& headers,
                             const std::vector<double>& values)
{
    if (headers.size() != values.size()) {
        throw std::runtime_error(
            "ResultWriter::write_row: headers and values size mismatch");
    }

    std::ofstream out(filename_, 
                      append_ ? std::ios::app : std::ios::trunc);
    if (!out.is_open()) {
        throw std::runtime_error(
            "ResultWriter::write_row: could not open output file");
    }

    // --------------------------------------------------
    // Write header if file is empty or not appending
    // --------------------------------------------------
    out.seekp(0, std::ios::end);
    bool empty = (out.tellp() == 0);

    if (empty) {
        // writer.write_comment("SSE simulation");
        // writer.write_comment(
        //     "Lx=" + std::to_string(prm.Lx) +
        //     " Ly=" + std::to_string(prm.Ly) +
        //     " beta=" + std::to_string(prm.beta)
        // );
        for (const auto& h : headers) {
            out << std::setw(16) << h;
        }
        out << '\n';
    }

    // --------------------------------------------------
    // Write values
    // --------------------------------------------------
    //out << std::scientific << std::setprecision(10);
    for (double v : values) {
        out << std::setw(16) << std::scientific << v;
    }
    out << '\n';
}

// --------------------------------------------------
// Write a comment line
// --------------------------------------------------
void ResultWriter::write_comment(const std::string& comment)
{
    std::ofstream out(filename_, std::ios::app);
    if (!out.is_open()) {
        throw std::runtime_error(
            "ResultWriter::write_comment: could not open output file");
    }
    out << "# " << comment << '\n';
}