from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        name="kmer_hashing.kmer_hashing",
        sources=["kmer_hashing/kmer_hashing.pyx"],
        extra_compile_args=['-O2', '-fPIC'],
    ),
]

setup(
    name="dedup_rate_predict",
    version="1.0",
    author="Wonjun Lim",
    author_email="cerutx@gmail.com",
    description="A package for Raw Data Quality Control",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/wjlim/dedup_rate_predict",
    packages=find_packages(),
    ext_modules=cythonize(extensions),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'dedup_rate_predict=dedup_rate_predict.dedup_rate_predict:main',
        ],
    },
    install_requires=[
        "Cython",
        "psutil",
    ],
)
