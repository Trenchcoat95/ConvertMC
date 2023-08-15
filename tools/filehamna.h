// Yea does that file even exist?  Probably doesn't work for xroot access.
inline bool filehamna(const std::string& filename) {
    struct stat buf;
    int retval = stat(filename.c_str(), &buf);
    if (retval==0 && !S_ISREG(buf.st_mode)) {
        cout << "filehamna(string) : Not a regular file" << endl;
    }
    // -1 is the error condition, meaning that the file isn't there
    return (retval == -1);
}
