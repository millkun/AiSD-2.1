import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import os

def bwt_transform(block):
    """Burrows-Wheeler Transform"""
    n = len(block)
    rotations = [block[i:] + block[:i] for i in range(n)]
    rotations.sort()
    transformed = bytes([rot[-1] for rot in rotations])
    index = rotations.index(block)
    return transformed, index

def mtf_transform(data):
    """Move-To-Front Transform"""
    alphabet = list(range(256))
    result = []
    for c in data:
        idx = alphabet.index(c)
        result.append(idx)
        # Move to front
        del alphabet[idx]
        alphabet.insert(0, c)
    return bytes(result)

def calculate_entropy(data):
    """Calculate Shannon entropy of data"""
    if not data:
        return 0
    counts = Counter(data)
    probs = [count/len(data) for count in counts.values()]
    return -sum(p * np.log2(p) for p in probs)

def process_file(filename, block_sizes):
    """Process file with different block sizes and calculate entropy"""
    file_size = os.path.getsize(filename)
    entropies = []
    
    with open(filename, 'rb') as f:
        for block_size in block_sizes:
            total_entropy = 0
            blocks_processed = 0
            f.seek(0)
            
            while True:
                block = f.read(block_size)
                if not block:
                    break
                
                # BWT
                bwt_data, _ = bwt_transform(block)
                
                # MTF
                mtf_data = mtf_transform(bwt_data)
                
                # Рассчитываю энтропию
                entropy = calculate_entropy(mtf_data)
                total_entropy += entropy * len(block)
                blocks_processed += len(block)
            
            if blocks_processed > 0:
                avg_entropy = total_entropy / blocks_processed
                entropies.append(avg_entropy)
                print(f"Block size: {block_size}, Entropy: {avg_entropy:.4f} bits/byte")
            else:
                entropies.append(0)
    
    return entropies

def plot_results(block_sizes, entropies):
    """График отношения"""
    plt.figure(figsize=(10, 6))
    plt.plot(block_sizes, entropies, 'o-')
    plt.xscale('log')
    plt.xticks(block_sizes, labels=[str(bs) for bs in block_sizes])
    plt.xlabel('Размер блоков')
    plt.ylabel('Энтропия')
    plt.title('Зависимость энтропии от размера блоков на enwik7.txt')
    plt.grid(True, which="both", ls="-")
    plt.tight_layout()
    plt.savefig('entropy_vs_block_size.png')
    plt.show()

filename = 'C:\\Users\\jeb32\\OneDrive\\Рабочий стол\\Всякое\\АиСД\\2 Курс; 2 Семестр\\Лабораторная1\\enwik7.txt'  # Make sure this file exists in the directory
block_sizes = [64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072]

entropies = process_file(filename, block_sizes)

plot_results(block_sizes, entropies)