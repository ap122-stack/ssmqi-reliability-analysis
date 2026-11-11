function digest = ssmqi_sha256_file(filePath)
% Return the lowercase SHA-256 digest of a file without modifying it.

if ~isfile(filePath)
    error('Cannot hash missing file: %s',filePath);
end
fid = fopen(filePath,'rb');
if fid < 0
    error('Cannot open file for hashing: %s',filePath);
end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid,Inf,'*uint8');
digest = sha256_bytes(bytes);
end

function digest = sha256_bytes(bytes)
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(bytes);
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end
