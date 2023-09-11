#include <crypto++/cryptlib.h>
#include <crypto++/files.h>
#include <crypto++/filters.h>
#include <crypto++/hex.h>
#include <crypto++/osrng.h>
#include <crypto++/sha.h>
#include <crypto++/eccrypto.h>
#include <crypto++/oids.h>
#include <crypto++/secblock.h>

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
using namespace CryptoPP;


typedef ECDSA<ECP, SHA256>::PrivateKey ECDSA_privateKey;
typedef ECDSA<ECP, SHA256>::PublicKey ECDSA_publicKey;
typedef ECDSA<ECP, SHA256>::Signer ECDSA_Signer;
typedef ECDSA<ECP, SHA256>::Verifier ECDSA_Verifier;

typedef ECIES<ECP>::PrivateKey ECC_privateKey;
typedef ECIES<ECP>::PublicKey ECC_publicKey;
typedef ECIES<ECP>::Encryptor ECC_Encryptor;
typedef ECIES<ECP>::Decryptor ECC_Decryptor;

class DEV_ECDH {
public:
	void GenerateKeys(unsigned int uiKeySize, std::string& outPrivateKey, std::string& outPublicKey) {
		using namespace CryptoPP;
		AutoSeededRandomPool prng(false, 256);
		ECC_privateKey pk;
		ECC_publicKey pu;

		pk.Initialize(prng, ASN1::secp256k1());
		pk.MakePublicKey(pu);

		ECC_Encryptor encryptor(pu);

	    HexEncoder pubEncoder(new StringSink(outPublicKey));
	    pu.DEREncode(pubEncoder);
	    pubEncoder.MessageEnd();

	   	ECC_Decryptor decryptor(pk);

	    HexEncoder prvEncoder(new StringSink(outPrivateKey));
	    pk.DEREncode(prvEncoder);
	    prvEncoder.MessageEnd();
		
	}

	std::string EncryptMessage(const std::string& outPublicKey, const std::string& sMsgToEncrypt) {
		using namespace CryptoPP;
	    
	    StringSource pubString(outPublicKey, true, new HexDecoder);
	    ECC_Encryptor encryptor(pubString);

	    size_t uiCipherTextSize = encryptor.CiphertextLength(sMsgToEncrypt.size());
	    std::string sCipherText;
	    sCipherText.resize(uiCipherTextSize);
	    RandomPool rnd;
	    encryptor.Encrypt(rnd, (byte*)(sMsgToEncrypt.c_str()), sMsgToEncrypt.size(), (byte*)(sCipherText.data()));
	    return sCipherText;
	}

	std::string DecryptMessage(const std::string& outPrivateKey, const std::string& sMsgToDecrytp)
	{
	    using namespace CryptoPP;
	    StringSource privString(outPrivateKey, true, new HexDecoder);
	   	ECC_Decryptor decryptor(privString);
	   	// Aqui estava com o tipo 'auto', mas tive que trocar para 'int', para que funcione com a versão c++ do docker ubuntu-14.04
	    unsigned int sPlainTextLen = decryptor.MaxPlaintextLength(sMsgToDecrytp.size());
	    std::string sDecryText;
	    sDecryText.resize(sPlainTextLen);
	    RandomPool rnd;
	    decryptor.Decrypt(rnd, (byte*)sMsgToDecrytp.c_str(), sMsgToDecrytp.size(), (byte*)sDecryText.data());
	    return sDecryText;
	}


};

class DEV_ECDSA {
	
public:
	void GenerateKeys(std::string& outPrivateKey, std::string& outPublicKey) {
		// Cria a chave privada 
		AutoSeededRandomPool prng(false, 256);
		ECDSA_privateKey pk;
		ECDSA_publicKey pu;

		pk.Initialize(prng, ASN1::secp256k1());
		pk.MakePublicKey(pu);

		ECDSA_Signer signer(pk);
	    
	    HexEncoder prvEncoder(new StringSink(outPrivateKey));
	    pk.DEREncode(prvEncoder);
	    prvEncoder.MessageEnd();

	    ECDSA_Verifier verifier(pu);

	    HexEncoder pubEncoder(new StringSink(outPublicKey));
	    pu.DEREncode(pubEncoder);
	    pubEncoder.MessageEnd();

	}
	
	std::string SignMessage(const std::string& outPrivateKey, const std::string message) {
		StringSource privString(outPrivateKey, true, new HexDecoder);
		ECDSA_Signer signer(privString);

		size_t siglen = signer.MaxSignatureLength();
		std::string signature(siglen, 0x00);
		AutoSeededRandomPool prng;

		siglen = signer.SignMessage( prng, (const byte*)&message[0], message.size(), (byte*)&signature[0]);
    	signature.resize(siglen);
    	
    	return signature;
	}

	void VerifyMessage(const std::string& outPublicKey, const std::string message, std::string signature) {
		StringSource pubString(outPublicKey, true, new HexDecoder);
		ECDSA_Verifier verifier(pubString);

		bool result = verifier.VerifyMessage((const byte*)&message[0], message.size(), (const byte*)&(signature[0]), signature.size());
		if( !result ) {
        	std::cout << "Falha ao verificar assinatura da mensagem" << std::endl;
    	} else {
        	std::cout << "Mensagem autenticada com sucesso\n" << std::endl;
    	}
	}

};
