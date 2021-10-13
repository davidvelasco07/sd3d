template <typename T>
T* malloc_host(size_t N, T value=T()) {
    T* ptr = (T*)(malloc(N*sizeof(T)));
    std::fill(ptr, ptr+N, value);
    return ptr;
}

template <typename T>
void Write_arrays(T *array, int N, char const *name, T value=T()){
    std::ofstream output;
    output.open(name);
    output.write((char*)array, N*sizeof(T));
    output.close();
}

class comm_cpu{
    public:
        int rank_src;
        int rank_dst;

        double *buffer_send;
        int buffer_size;

        int build(int src, int dst, int size){
            rank_src = src;
            rank_dst = dst;
            buffer_send = malloc_host<double>(size);
            buffer_size = size;
            return 1;
        }	
};
