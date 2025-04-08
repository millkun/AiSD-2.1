import matplotlib.pyplot as plt
import os
from typing import List, Tuple, NamedTuple
from pathlib import Path
import time

class Token(NamedTuple):
    offset: int    # расстояние до начала найденного совпадения в окне
    length: int    # длина совпадения
    next_char: str # следующий символ после совпадения
    is_last: bool  # если True, то токен завершающий

def compress_lz77(data: str, window_size: int = 4096, lookahead_size: int = 16) -> List[Token]:
    tokens = []
    i = 0
    n = len(data)
    
    while i < n:
        best_length = 0
        best_offset = 0
        start = max(0, i - window_size)
        
        # Поиск наилучшего совпадения в окне
        for j in range(start, i):
            length = 0
            while (length < lookahead_size and 
                   i + length < n and 
                   j + length < i and 
                   data[j + length] == data[i + length]):
                length += 1
            
            if length > best_length:
                best_length = length
                best_offset = i - j
        
        # Создание токена
        if i + best_length < n:
            next_char = data[i + best_length]
            token = Token(best_offset, best_length, next_char, False)
            tokens.append(token)
            i += best_length + 1
        else:
            token = Token(best_offset, best_length, '', True)
            tokens.append(token)
            i += best_length
    
    return tokens

def calculate_compression_ratio(original_size: int, tokens: List[Token]) -> float:
    # Размер токена: offset (2 байта), length (1 байт), next_char (1 байт), is_last (1 байт)
    token_size = 5
    compressed_size = len(tokens) * token_size
    return original_size / compressed_size if compressed_size > 0 else 1

def analyze_buffer_impact(data: str, max_window_size: int = 32768, step: int = 1024) -> Tuple[List[int], List[float], List[float]]:
    window_sizes = []
    ratios = []
    times = []
    
    for window_size in range(1024, max_window_size + 1, step):
        start_time = time.time()
        tokens = compress_lz77(data, window_size=window_size)
        elapsed = time.time() - start_time
        
        ratio = calculate_compression_ratio(len(data), tokens)
        window_sizes.append(window_size)
        ratios.append(ratio)
        times.append(elapsed)
        
        print(f"Window: {window_size//1024}KB, Ratio: {ratio:.2f}, Time: {elapsed:.2f}s")
    
    return window_sizes, ratios, times

def plot_results(window_sizes: List[int], ratios: List[float], times: List[float]):
    plt.figure(figsize=(12, 8))
    
    # График коэффициента сжатия
    plt.subplot(2, 1, 1)
    plt.plot(window_sizes, ratios, marker='o', linestyle='-', color='b')
    plt.title('Зависимость коэффициента сжатия от размера буфера (LZ77 на enwik7)')
    plt.xlabel('Размер буфера (окна) в байтах')
    plt.ylabel('Коэффициент сжатия')
    plt.grid(True)
    
    # График времени выполнения
    plt.subplot(2, 1, 2)
    plt.plot(window_sizes, times, marker='o', linestyle='-', color='r')
    plt.title('Зависимость времени сжатия от размера буфера')
    plt.xlabel('Размер буфера (окна) в байтах')
    plt.ylabel('Время сжатия (секунды)')
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

def read_enwik7(file_path: str) -> str:
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        return f.read()

if __name__ == "__main__":
    # Укажите путь к файлу enwik7
    enwik7_path = 'C:\\Users\\jeb32\\OneDrive\\Рабочий стол\\Всякое\\АиСД\\2 Курс; 2 Семестр\\Лабораторная1\\enwik7.txt'
    
    if not Path(enwik7_path).exists():
        print(f"Файл {enwik7_path} не найден!")
        print("Скачайте его с http://mattmahoney.net/dc/textdata.html")
        exit(1)
    
    print("Чтение файла enwik7...")
    data = read_enwik7(enwik7_path)
    print(f"Размер данных: {len(data)/1024:.1f} KB")
    
    # Для теста можно взять часть данных (полные данные требуют много времени)
    sample_size = min(100000, len(data))  # Можно увеличить для более точных результатов
    test_data = data[:sample_size]
    
    print("\nАнализ влияния размера буфера...")
    window_sizes, ratios, times = analyze_buffer_impact(
        test_data, 
        max_window_size=16384,  # Максимальный размер буфера для анализа
        step=1024               # Шаг изменения размера буфера
    )
    
    print("\nПостроение графиков...")
    plot_results(window_sizes, ratios, times)