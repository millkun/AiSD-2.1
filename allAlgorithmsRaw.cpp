#include <iostream>
#include <string>
#include <queue>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <tuple>
#include <sstream>
#include <cstring>
#include <memory>
#include <functional>
#include <map>

using namespace std;

// ----------------------------------- HA [Huffman Algorithm] -----------------------------------

struct TreeNode {
    char character;
    unsigned frequency, code;
    TreeNode* left, * right;

    bool operator>(const TreeNode& b) const {
        return frequency > b.frequency;
    }

    TreeNode(char chara, int freq)
        : character(chara), frequency(freq), code(0u), left(nullptr), right(nullptr) {
    }
    TreeNode() : character(0), frequency(0), code(0u), left(nullptr), right(nullptr) {}
};

TreeNode* buildHuffmanTree(const map<char, unsigned>& frequencies) {
    auto cmp = [](const TreeNode* a, const TreeNode* b) { return *a > *b; };
    priority_queue<TreeNode*, vector<TreeNode*>, decltype(cmp)> pairs(cmp);

    for (const auto& pair : frequencies) {
        pairs.push(new TreeNode(pair.first, pair.second));
    }

    while (pairs.size() > 1) {
        TreeNode* pair1 = pairs.top();
        pairs.pop();
        TreeNode* pair2 = pairs.top();
        pairs.pop();

        TreeNode* temp = new TreeNode('*', pair1->frequency + pair2->frequency);
        temp->left = pair1;
        temp->right = pair2;
        pairs.push(temp);
    }

    return pairs.empty() ? nullptr : pairs.top();
}

void assignCodes(map<char, string>& codes, TreeNode* node, string code) {
    if (!node->left && !node->right) {
        codes[node->character] = code;
        return;
    }
    if (node->left) assignCodes(codes, node->left, code + "0");
    if (node->right) assignCodes(codes, node->right, code + "1");
}

void writeTree(ofstream& out, TreeNode* node) {
    if (!node->left && !node->right) {
        out.put(1);
        out.put(node->character);
        return;
    }
    out.put(0);
    if (node->left) writeTree(out, node->left);
    if (node->right) writeTree(out, node->right);
}

TreeNode* readTree(ifstream& in) {
    char flag = in.get();
    if (flag == 1) {
        char ch = in.get();
        return new TreeNode(ch, 0);
    }
    TreeNode* node = new TreeNode('*', 0);
    node->left = readTree(in);
    node->right = readTree(in);
    return node;
}

void compressHA(const string& inputFile) {
    ifstream in(inputFile, ios::binary);
    if (!in) {
        cerr << "Error opening input file\n";
        return;
    }

    // First pass: calculate frequencies
    map<char, unsigned> frequencies;
    char ch;
    while (in.get(ch)) {
        frequencies[ch]++;
    }
    in.clear();
    in.seekg(0);

    TreeNode* root = buildHuffmanTree(frequencies);
    map<char, string> codes;
    assignCodes(codes, root, "");

    ofstream out("compressed-HA.bin", ios::binary);
    if (!out) {
        cerr << "Error creating output file\n";
        return;
    }

    // Write Huffman tree
    writeTree(out, root);

    // Write compressed data
    string buffer;
    while (in.get(ch)) {
        buffer += codes[ch];
        while (buffer.size() >= 8) {
            bitset<8> bits(buffer.substr(0, 8));
            out.put(static_cast<char>(bits.to_ulong()));
            buffer = buffer.substr(8);
        }
    }

    // Write remaining bits
    if (!buffer.empty()) {
        while (buffer.size() < 8) buffer += '0';
        bitset<8> bits(buffer);
        out.put(static_cast<char>(bits.to_ulong()));
    }

    in.close();
    out.close();
}

void decompressHA(const string& compressedFile) {
    ifstream in(compressedFile, ios::binary);
    if (!in) {
        cerr << "Error opening compressed file\n";
        return;
    }

    TreeNode* root = readTree(in);

    ofstream out("decompressed-HA.txt", ios::binary);
    if (!out) {
        cerr << "Error creating output file\n";
        return;
    }

    TreeNode* current = root;
    char byte;
    while (in.get(byte)) {
        bitset<8> bits(byte);
        for (int i = 7; i >= 0; --i) {
            if (bits[i]) {
                current = current->right;
            }
            else {
                current = current->left;
            }

            if (!current->left && !current->right) {
                out.put(current->character);
                current = root;
            }
        }
    }

    in.close();
    out.close();
}

// ----------------------------------- RLE [Run-Length Encoding] -----------------------------------

vector<char> readFile(const string& filename) {
    ifstream file(filename, ios::binary);
    if (!file.is_open()) throw runtime_error("Cannot open file: " + filename);

    file.seekg(0, ios::end);
    size_t fileSize = file.tellg();
    file.seekg(0, ios::beg);

    vector<char> buffer(fileSize);
    if (fileSize > 0) file.read(buffer.data(), fileSize);
    return buffer;
}

void writeFile(const string& filename, const vector<char>& data) {
    ofstream file(filename, ios::binary);
    if (!file.is_open()) throw runtime_error("Cannot create file: " + filename);
    if (!data.empty()) file.write(data.data(), data.size());
}

vector<char> compressRLE(const vector<char>& input) {
    vector<char> compressed;
    if (input.empty()) return compressed;

    size_t i = 0;
    const size_t n = input.size();

    while (i < n) {
        char current = input[i];
        size_t count = 1;

        // Подсчет повторяющихся символов
        while (i + count < n && input[i + count] == current && count < 127) {
            count++;
        }

        if (count > 1) {
            compressed.push_back(static_cast<char>(count));
            compressed.push_back(current);
            i += count;
        }
        else {
            // Подсчет уникальных символов
            size_t unique_count = 0;
            while (i + unique_count < n &&
                (i + unique_count + 1 >= n || input[i + unique_count] != input[i + unique_count + 1]) &&
                unique_count < 127) {
                unique_count++;
            }

            if (unique_count == 0) unique_count = 1; // Защита от бесконечного цикла

            compressed.push_back(static_cast<char>(-static_cast<signed char>(unique_count)));
            for (size_t j = 0; j < unique_count && i + j < n; j++) {
                compressed.push_back(input[i + j]);
            }
            i += unique_count;
        }
    }

    return compressed;
}

vector<char> decompressRLE(const vector<char>& compressed) {
    vector<char> decompressed;
    if (compressed.empty()) return decompressed;

    size_t i = 0;
    const size_t n = compressed.size();

    while (i < n) {
        signed char count = static_cast<signed char>(compressed[i]);

        if (count > 0) {
            // Обработка повторяющихся символов
            if (i + 1 >= n) break; // Недостаточно данных
            char value = compressed[i + 1];
            for (int j = 0; j < count; j++) {
                decompressed.push_back(value);
            }
            i += 2;
        }
        else {
            // Обработка уникальных символов
            int unique_count = -count;
            if (i + unique_count >= n) break; // Недостаточно данных

            for (int j = 1; j <= unique_count && i + j < n; j++) {
                decompressed.push_back(compressed[i + j]);
            }
            i += 1 + unique_count;
        }
    }

    return decompressed;
}

void compressRLE(const string& inputFile) {
    try {
        cout << "Reading file: " << inputFile << endl;
        vector<char> inputData = readFile(inputFile);
        cout << "File size: " << inputData.size() << " bytes" << endl;

        vector<char> compressedData = compressRLE(inputData);
        writeFile("compressed-RLE.bin", compressedData);

        cout << "RLE compression successful. Compressed size: "
            << compressedData.size() << " bytes (ratio: "
            << (inputData.empty() ? 0 : (double)compressedData.size() / inputData.size() * 100)
            << "%)" << endl;
    }
    catch (const exception& e) {
        cerr << "Error during RLE compression: " << e.what() << endl;
    }
}

void decompressRLE(const string& inputFile) {
    try {
        cout << "Reading compressed file: " << inputFile << endl;
        vector<char> compressedData = readFile(inputFile);
        cout << "Compressed size: " << compressedData.size() << " bytes" << endl;

        vector<char> decompressedData = decompressRLE(compressedData);
        writeFile("decompressed-RLE.txt", decompressedData);

        cout << "RLE decompression successful. Decompressed size: "
            << decompressedData.size() << " bytes" << endl;
    }
    catch (const exception& e) {
        cerr << "Error during RLE decompression: " << e.what() << endl;
    }
}

// ----------------------------------- BWT [Burrows–Wheeler Transform] -----------------------------------

// Размер блока для обработки (можно настроить под вашу систему)
const size_t BLOCK_SIZE = 4096; // 8KB | 16384
const char EOF_MARKER = 0x01;

// BWT functions
void bwtencode(const string& inputFile, const string& outputFile);
void bwtdecode(const string& inputFile, const string& outputFile);
int C(char c, const string& L);
int Occ(char c, int r, const string& L);

void compressBWT(const string& inputFile) {
    string outputFile = "compressed-BWT.bin";
    bwtencode(inputFile, outputFile);
    cout << "BWT encoding completed. Result saved to " << outputFile << endl;
}

void decompressBWT(const string& inputFile) {
    string outputFile = "decompressed-BWT.txt";
    bwtdecode(inputFile, outputFile);
    cout << "BWT decoding completed. Result saved to " << outputFile << endl;
}

void bwtencode(const string& inputFile, const string& outputFile) {
    ifstream inFile(inputFile, ios::binary);
    ofstream outFile(outputFile, ios::binary);

    if (!inFile) {
        cerr << "Error opening input file!" << endl;
        return;
    }
    if (!outFile) {
        cerr << "Error opening output file!" << endl;
        return;
    }

    char EOF_MARKER = 0x01;
    size_t total_processed = 0;
    size_t block_num = 0;

    while (true) {
        // Читаем блок данных
        vector<char> buffer(BLOCK_SIZE);
        inFile.read(buffer.data(), BLOCK_SIZE);
        size_t bytes_read = inFile.gcount();

        if (bytes_read == 0 && block_num > 0) break;
        if (bytes_read == 0) {
            cerr << "Input file is empty!" << endl;
            return;
        }

        // Добавляем маркер конца блока (разные для каждого блока)
        char block_marker = EOF_MARKER + block_num;
        buffer.resize(bytes_read);
        buffer.push_back(block_marker);

        // Создаем rotations для текущего блока
        string text(buffer.begin(), buffer.end());
        string L;
        string last_str = text;
        size_t len_text = text.length();
        vector<string> rotate_list;
        rotate_list.push_back(last_str);

        for (size_t i = 0; i < len_text - 1; ++i) {
            last_str = last_str.back() + last_str.substr(0, len_text - 1);
            rotate_list.push_back(last_str);
        }

        // Сортируем и строим BWT
        sort(rotate_list.begin(), rotate_list.end());
        for (const auto& eachstr : rotate_list) {
            L += eachstr.back();
        }

        // Записываем размер блока и сам блок
        uint32_t block_size = L.size();
        outFile.write(reinterpret_cast<const char*>(&block_size), sizeof(block_size));
        outFile.write(L.c_str(), L.size());

        total_processed += bytes_read;
        block_num++;
        cout << "Processed block " << block_num << " (" << bytes_read << " bytes)" << endl;
    }

    cout << "Encoded total " << total_processed << " bytes in " << block_num << " blocks." << endl;
}

void bwtdecode(const string& inputFile, const string& outputFile) {
    ifstream inFile(inputFile, ios::binary);
    ofstream outFile(outputFile, ios::binary);

    if (!inFile) {
        cerr << "Error opening input file!" << endl;
        return;
    }
    if (!outFile) {
        cerr << "Error opening output file!" << endl;
        return;
    }

    size_t block_num = 0;
    size_t total_processed = 0;

    while (true) {
        // Читаем размер блока
        uint32_t block_size;
        inFile.read(reinterpret_cast<char*>(&block_size), sizeof(block_size));

        if (inFile.gcount() == 0) break; // Конец файла

        // Читаем сам блок
        vector<char> L_buffer(block_size);
        inFile.read(L_buffer.data(), block_size);

        if (inFile.gcount() != block_size) {
            cerr << "Error reading block " << block_num << endl;
            return;
        }

        string L(L_buffer.begin(), L_buffer.end());

        // Определяем маркер конца блока
        char block_marker = EOF_MARKER + block_num;
        size_t marker_pos = L.find(block_marker);
        if (marker_pos == string::npos) {
            cerr << "Error: EOF marker not found in block " << block_num << endl;
            return;
        }

        // Восстанавливаем блок
        vector<int> counts(256, 0);
        for (char c : L) {
            counts[static_cast<unsigned char>(c)]++;
        }

        vector<int> C_table(256, 0);
        for (int i = 1; i < 256; ++i) {
            C_table[i] = C_table[i - 1] + counts[i - 1];
        }

        vector<int> LF(L.size());
        vector<int> next_pos(256, 0);
        for (size_t i = 0; i < L.size(); ++i) {
            unsigned char c = L[i];
            LF[i] = C_table[c] + next_pos[c];
            next_pos[c]++;
        }

        string text;
        int row = marker_pos;
        for (size_t i = 0; i < L.size() - 1; ++i) {
            text.push_back(L[row]);
            row = LF[row];
        }

        reverse(text.begin(), text.end());
        outFile.write(text.c_str(), text.size());
        total_processed += text.size();
        block_num++;

        cout << "Decoded block " << block_num << " (" << text.size() << " bytes)" << endl;
    }

    cout << "Decoded total " << total_processed << " bytes from " << block_num << " blocks." << endl;
}

int C(char c, const string& L) {
    // Count all characters in L whose lexicographic order is less than "c"
    int num_Tc = 0;
    for (char eachchr : L) {
        if (c > eachchr) {
            num_Tc++;
        }
    }
    return num_Tc - 1;
}

int Occ(char c, int r, const string& L) {
    // Count the number of characters c before line r in L
    int row = -1;
    int num_Lc = 0;
    if (r == 0) {
        return 0;
    }
    for (char eachchr : L) {
        row++;
        if (eachchr == c) {
            num_Lc++;
        }
        if (row == r - 1) {
            break;
        }
    }
    return num_Lc;
}

// ----------------------------------- LZ77 [Lempel-Ziv 1977] -----------------------------------

// Структура для хранения закодированных данных LZ77
struct Node {
    int offset;
    int length;
    char next;
};

// Функция для поиска совпадения в буфере
pair<int, int> findMatching(const vector<char>& searchBuffer, const vector<char>& lookAheadBuffer) {
    int maxLength = 0;
    int bestOffset = 0;

    // Ищем максимальное совпадение между searchBuffer и lookAheadBuffer
    for (int i = 0; i < searchBuffer.size(); ++i) {
        int currentLength = 0;
        while (currentLength < lookAheadBuffer.size() &&
            i + currentLength < searchBuffer.size() &&
            searchBuffer[i + currentLength] == lookAheadBuffer[currentLength]) {
            currentLength++;
        }

        if (currentLength > maxLength) {
            maxLength = currentLength;
            bestOffset = searchBuffer.size() - i;
        }
    }

    return make_pair(bestOffset, maxLength);
}

// Функция сжатия LZ77
void compressLZ77(const string& inputFile) {
    ifstream inFile(inputFile, ios::binary);
    ofstream outFile("compressed-LZ77.bin", ios::binary);

    if (!inFile || !outFile) {
        cerr << "Error opening files!" << endl;
        return;
    }

    const size_t searchBufferSize = 1024; // Размер буфера поиска
    const size_t lookAheadBufferSize = 256; // Размер буфера просмотра вперед

    vector<char> searchBuffer;
    vector<char> lookAheadBuffer;
    list<Node> encoded;

    // Читаем первые данные в буфер просмотра вперед
    lookAheadBuffer.resize(lookAheadBufferSize);
    inFile.read(lookAheadBuffer.data(), lookAheadBufferSize);
    size_t bytesRead = inFile.gcount();
    lookAheadBuffer.resize(bytesRead);

    while (!lookAheadBuffer.empty()) {
        // Ищем совпадение
        auto match = findMatching(searchBuffer, lookAheadBuffer);
        int offset = match.first;
        int length = match.second;

        // Читаем следующий символ
        char nextChar = (length < lookAheadBuffer.size()) ? lookAheadBuffer[length] : '\0';

        // Добавляем в закодированные данные
        encoded.push_back({ offset, length, nextChar });

        // Обновляем буферы
        size_t shiftAmount = length + 1;

        // Добавляем обработанные данные в searchBuffer
        for (size_t i = 0; i < shiftAmount && i < lookAheadBuffer.size(); ++i) {
            searchBuffer.push_back(lookAheadBuffer[i]);
        }

        // Удаляем лишние данные из searchBuffer, если он превысил размер
        if (searchBuffer.size() > searchBufferSize) {
            size_t excess = searchBuffer.size() - searchBufferSize;
            searchBuffer.erase(searchBuffer.begin(), searchBuffer.begin() + excess);
        }

        // Удаляем обработанные данные из lookAheadBuffer
        lookAheadBuffer.erase(lookAheadBuffer.begin(), lookAheadBuffer.begin() + min(shiftAmount, lookAheadBuffer.size()));

        // Если lookAheadBuffer пуст, читаем новые данные
        if (lookAheadBuffer.empty() || lookAheadBuffer.size() < lookAheadBufferSize) {
            vector<char> newData(lookAheadBufferSize - lookAheadBuffer.size());
            inFile.read(newData.data(), newData.size());
            size_t newBytes = inFile.gcount();
            newData.resize(newBytes);

            lookAheadBuffer.insert(lookAheadBuffer.end(), newData.begin(), newData.end());
        }
    }

    // Записываем закодированные данные в файл
    for (const auto& node : encoded) {
        outFile.write(reinterpret_cast<const char*>(&node.offset), sizeof(node.offset));
        outFile.write(reinterpret_cast<const char*>(&node.length), sizeof(node.length));
        outFile.write(&node.next, sizeof(node.next));
    }

    inFile.close();
    outFile.close();
}

// Функция распаковки LZ77
void decompressLZ77(const string& inputFile) {
    ifstream inFile(inputFile, ios::binary);
    ofstream outFile("decompressed-LZ77.txt", ios::binary);

    if (!inFile || !outFile) {
        cerr << "Error opening files!" << endl;
        return;
    }

    vector<char> decodedData;
    list<Node> encoded;

    // Читаем закодированные данные из файла
    while (true) {
        Node node;
        inFile.read(reinterpret_cast<char*>(&node.offset), sizeof(node.offset));
        inFile.read(reinterpret_cast<char*>(&node.length), sizeof(node.length));
        inFile.read(&node.next, sizeof(node.next));

        if (inFile.eof()) break;

        encoded.push_back(node);
    }

    // Декодируем данные
    for (const auto& node : encoded) {
        if (node.length > 0) {
            // Копируем последовательность из уже декодированных данных
            int start = decodedData.size() - node.offset;
            for (int i = 0; i < node.length; ++i) {
                if (start + i >= 0 && start + i < decodedData.size()) {
                    decodedData.push_back(decodedData[start + i]);
                }
            }
        }

        // Добавляем следующий символ
        if (node.next != '\0') {
            decodedData.push_back(node.next);
        }
    }

    // Записываем декодированные данные в файл
    outFile.write(decodedData.data(), decodedData.size());

    inFile.close();
    outFile.close();
}

// ----------------------------------- LZW [Lempel-Ziv-Welch] -----------------------------------

const int MAX_DICT_SIZE = 4096; // Максимальный размер словаря

// Функция для сжатия данных с использованием алгоритма LZW
void compressLZW(const string& inputFile) {
    ifstream inFile(inputFile, ios::binary);
    ofstream outFile("compressed-LZW.bin", ios::binary);

    if (!inFile.is_open() || !outFile.is_open()) {
        cerr << "Error opening files!" << endl;
        return;
    }

    unordered_map<string, int> dictionary;
    for (int i = 0; i < 256; i++) {
        dictionary[string(1, char(i))] = i;
    }

    string current = "";
    int dictSize = 256;
    vector<int> output;

    char byte;
    while (inFile.get(byte)) {
        string next = current + byte;
        if (dictionary.count(next)) {
            current = next;
        }
        else {
            output.push_back(dictionary[current]);
            if (dictSize < MAX_DICT_SIZE) {
                dictionary[next] = dictSize++;
            }
            current = string(1, byte);
        }
    }

    if (!current.empty()) {
        output.push_back(dictionary[current]);
    }

    // Запись сжатых данных в файл
    for (int code : output) {
        outFile.write(reinterpret_cast<const char*>(&code), sizeof(code));
    }

    inFile.close();
    outFile.close();
}

// Функция для разжатия данных с использованием алгоритма LZW
void decompressLZW(const string& inputFile) {
    ifstream inFile(inputFile, ios::binary);
    ofstream outFile("decompressed-LZW.txt", ios::binary);

    if (!inFile.is_open() || !outFile.is_open()) {
        cerr << "Error opening files!" << endl;
        return;
    }

    unordered_map<int, string> dictionary;
    for (int i = 0; i < 256; i++) {
        dictionary[i] = string(1, char(i));
    }

    int dictSize = 256;
    int previousCode;
    inFile.read(reinterpret_cast<char*>(&previousCode), sizeof(previousCode));
    string current = dictionary[previousCode];
    outFile << current;

    int code;
    while (inFile.read(reinterpret_cast<char*>(&code), sizeof(code))) {
        string entry;
        if (dictionary.count(code)) {
            entry = dictionary[code];
        }
        else if (code == dictSize) {
            entry = current + current[0];
        }
        else {
            cerr << "Invalid LZW code: " << code << endl;
            break;
        }

        outFile << entry;

        // Добавление новой строки в словарь
        if (dictSize < MAX_DICT_SIZE) {
            dictionary[dictSize++] = current + entry[0];
        }

        current = entry;
    }

    inFile.close();
    outFile.close();
}

// ----------------------------------- MTF [Move-To-Front] -----------------------------------

vector<unsigned char> symbolTable(256);

void moveToFront(int i) {
    unsigned char t = symbolTable[i];
    for (int z = i - 1; z >= 0; z--)
        symbolTable[z + 1] = symbolTable[z];
    symbolTable[0] = t;
}

void fillSymbolTable() {
    for (int x = 0; x < 256; x++)
        symbolTable[x] = static_cast<unsigned char>(x);
}

string mtfEncode(const string& str) {
    fillSymbolTable();
    vector<int> output;
    for (unsigned char c : str) {
        for (int i = 0; i < 256; i++) {
            if (c == symbolTable[i]) {
                output.push_back(i);
                moveToFront(i);
                break;
            }
        }
    }
    string r;
    for (int num : output) {
        r += to_string(num) + " ";
    }
    return r;
}

string mtfDecode(const string& str) {
    fillSymbolTable();
    istringstream iss(str);
    vector<int> output;
    copy(istream_iterator<int>(iss), istream_iterator<int>(), back_inserter(output));

    string r;
    for (int num : output) {
        r += static_cast<char>(symbolTable[num]);
        moveToFront(num);
    }
    return r;
}

void compressMTF(const string& inputFile) {
    const size_t BUFFER_SIZE = 4096; // 4KB buffer
    ifstream inFile(inputFile, ios::binary);
    ofstream outFile("compressed-MTF.bin", ios::binary);

    if (!inFile) {
        cerr << "Error opening input file!" << endl;
        return;
    }

    if (!outFile) {
        cerr << "Error creating output file!" << endl;
        return;
    }

    vector<char> buffer(BUFFER_SIZE);
    while (inFile.read(buffer.data(), BUFFER_SIZE)) {
        string chunk(buffer.begin(), buffer.begin() + inFile.gcount());
        string encoded = mtfEncode(chunk);
        outFile.write(encoded.c_str(), encoded.size());
    }

    // Process remaining bytes if any
    if (inFile.gcount() > 0) {
        string chunk(buffer.begin(), buffer.begin() + inFile.gcount());
        string encoded = mtfEncode(chunk);
        outFile.write(encoded.c_str(), encoded.size());
    }

    cout << "MTF compression completed. Result saved to compressed-MTF.bin" << endl;
}

void decompressMTF(const string& inputFile) {
    const size_t BUFFER_SIZE = 4096; // 4KB buffer
    ifstream inFile(inputFile, ios::binary);
    ofstream outFile("decompressed-MTF.txt", ios::binary); // Изменил расширение на .bin для бинарных данных

    if (!inFile) {
        cerr << "Error opening input file!" << endl;
        return;
    }

    if (!outFile) {
        cerr << "Error creating output file!" << endl;
        return;
    }

    // Читаем весь файл (так как MTF коды - это текст)
    string encodedData;
    vector<char> buffer(BUFFER_SIZE);
    while (inFile.read(buffer.data(), BUFFER_SIZE)) {
        encodedData.append(buffer.begin(), buffer.begin() + inFile.gcount());
    }
    // Process remaining bytes if any
    if (inFile.gcount() > 0) {
        encodedData.append(buffer.begin(), buffer.begin() + inFile.gcount());
    }

    string decoded = mtfDecode(encodedData);
    outFile.write(decoded.c_str(), decoded.size());

    cout << "MTF decompression completed. Result saved to decompressed-MTF.bin" << endl;
}

int main() {
    int choice;
    string inputFile;

    std::cout << "Select compression algorithm:\n";
    std::cout << "1. HA - Done\n";
    std::cout << "2. BWT - Done\n";
    std::cout << "3. MTF - Done\n";
    std::cout << "4. RLE - Done\n";
    std::cout << "5. LZ77 - Done\n";
    std::cout << "6. LZW - Done\n";
    std::cout << "Enter your choice: ";
    std::cin >> choice;

    std::cout << "Enter input file name: ";
    std::cin >> inputFile;

    switch (choice) {
    case 1: // Correct
        compressHA(inputFile);
        decompressHA("compressed-HA.bin");
        break;
    case 2: // Correct
        compressBWT(inputFile);
        decompressBWT("compressed-BWT.bin");
        break;
    case 3: // Correct
        compressMTF(inputFile);
        decompressMTF("compressed-MTF.bin");
        break;
    case 4: // Correct
        compressRLE(inputFile);
        decompressRLE("compressed-RLE.bin");
        break;
    case 5: // Correct
        compressLZ77(inputFile);
        decompressLZ77("compressed-LZ77.bin");
        break;
    case 6: // Correct
        compressLZW(inputFile);
        decompressLZW("compressed-LZW.bin");
        break;
    default:
        cout << "Invalid choice!\n";
    }

    return 0;
}
