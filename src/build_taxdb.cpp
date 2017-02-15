/*
 * Copyright 2017, Florian Breitwieser
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "taxdb.h"

#include <iostream>
using namespace std;

int main(int argc, char **argv) {
    std::string database_dir = argv[0];
    TaxonomyDB taxdb;
    taxdb.writeTaxonomyIndex(
            std::cout,
            database_dir + "/taxonomy/nodes.dmp",
            database_dir + "/taxonomy/names.dmp");

}
