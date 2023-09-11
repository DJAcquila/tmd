#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include "ns3/address.h"
#include "util_consensus.h"
using namespace ns3;

class Block {
public:
	Block();
	Block(Transaction transaction, int decision, 
		int blockHeight, int contentID, 
		int parentBlockContentID, std::string parentBlockHash, 
		int byteblockSize, double timeCreated, double timeReceived, 
		Ipv4Address receivedFrom = Ipv4Address("0.0.0.0"), int receivedFromId = -1, 
		Ipv4Address sender = Ipv4Address("0.0.0.0"), int senderId = -1);
	Block(const Block &blockSrc);
	virtual ~Block();

	Transaction GetTransaction(void);
	void SetTransaction(Transaction t);
	
	int GetBlockHeight(void) const;
	void SetBlockHeight(int height);

	int GetContentID(void) const;
	void SetContentID(int contentID);

	int GetDecision(void) const;
	void SetDecision(int decision);

	int GetParentBlockContentID(void) const;
	void SetParentBlockContentID(int parentBlockContentID);

	std::string GetParentBlockHash(void) const;
	void SetParentBlockHash(std::string parentBlockHash);

	int GetByteblockSize(void) const;
	void SetByteblockSize(int byteblockSize);

	std::string GetBlockHash (void) const;
	std::string GenerateBLockHash(void);
	bool VerifyHash(std::string recv_hash);

	double GetTimeCreated(void) const;

	double GetTimeReceived(void) const;

	Ipv4Address GetReceivedFrom(void) const;
	void SetReceivedFrom(Ipv4Address receivedFrom);

	int GetReceivedFromId(void) const;
	void SetReceivedFromId(int receivedFromId);

	Ipv4Address GetSender(void) const;
	void SetSender(Ipv4Address sender);

	int GetSenderId(void) const;
	void SetSenderId(int senderId);

	bool IsParent(const Block &block) const;

	bool IsChild(const Block &block) const;

	void Print(void);

	Block& operator= (const Block &blockSrc) {
		m_blockHeight = blockSrc.m_blockHeight;
		m_contentID = blockSrc.m_contentID;
		m_parentBlockContentID = blockSrc.m_parentBlockContentID;
		m_byteblockSize = blockSrc.m_byteblockSize;
		m_timeCreated = blockSrc.m_timeCreated;
		m_timeReceived = blockSrc.m_timeReceived;
		m_receivedFrom = blockSrc.m_receivedFrom;
		m_transaction = blockSrc.m_transaction;
		m_decision = blockSrc.m_decision;
		m_receivedFromId = blockSrc.m_receivedFromId;
		m_parentBlockHash = blockSrc.m_parentBlockHash;
		m_senderId = blockSrc.m_senderId;
		m_sender = blockSrc.m_sender;

		return *this;
	}

	friend bool operator== (const Block &b1, const Block &b2);

	friend std::ostream& operator << (std::ostream &out, const Block &block);

protected:
	Transaction 	m_transaction;
	int 			m_blockHeight;
	int 			m_contentID;
	int 			m_byteblockSize;
	int 			m_parentBlockContentID;
	double 			m_timeCreated;
	double 			m_timeReceived;
	std::string 	m_blockHash;
	Ipv4Address 	m_receivedFrom;
	int 			m_receivedFromId;
	Ipv4Address 	m_sender;
	int 			m_senderId;
	int 			m_decision;
	std::string 	m_parentBlockHash;
};

Block::Block(Transaction transaction, int decision, int blockHeight, int contentID, 
	int parentBlockContentID, std::string parentBlockHash, 
	int byteblockSize, double timeCreated, double timeReceived, 
	Ipv4Address receivedFrom, int receivedFromId,
	Ipv4Address sender, int senderId){

	m_blockHeight = blockHeight;
	m_contentID = contentID;
	m_parentBlockContentID = parentBlockContentID;
	m_byteblockSize = byteblockSize;
	m_timeCreated = timeCreated;
	m_timeReceived = timeReceived;
	m_receivedFrom = receivedFrom;
	m_transaction = transaction;
	m_receivedFromId = receivedFromId;
	m_senderId = senderId;
	m_sender = sender;
	m_decision = decision;
	m_blockHash = GenerateBLockHash();
	m_parentBlockHash = parentBlockHash;
}

Block::Block () {
	Transaction t;
	Block(t, -1, -1, -1, -1, std::string(""), -1, -1.0, -1.0, Ipv4Address("0.0.0.0"), -1, Ipv4Address("0.0.0.0"), -1);
}

Block::Block (const Block& blockSrc) {
	m_blockHeight = blockSrc.m_blockHeight;
	m_contentID = blockSrc.m_contentID;
	m_parentBlockContentID = blockSrc.m_parentBlockContentID;
	m_byteblockSize = blockSrc.m_byteblockSize;
	m_timeCreated = blockSrc.m_timeCreated;
	m_timeReceived = blockSrc.m_timeReceived;
	m_receivedFrom = blockSrc.m_receivedFrom;
	m_receivedFromId = blockSrc.m_receivedFromId;
	m_transaction = blockSrc.m_transaction;
	m_decision = blockSrc.m_decision;
	m_parentBlockHash = blockSrc.m_parentBlockHash;
	m_senderId = blockSrc.m_senderId;
	m_sender = blockSrc.m_sender;
}

Block::~Block() {}

Transaction Block::GetTransaction(void){
	return m_transaction;
}

std::string Block::GetBlockHash(void) const {
	return m_blockHash;
}

std::string Block::GenerateBLockHash(void)  {
	std::string msg = m_transaction.getTxHash()
					+ std::to_string(m_decision)
					+ std::to_string(m_contentID)
					+ std::to_string(m_parentBlockContentID)
					+ m_parentBlockHash
					+ std::to_string(m_byteblockSize)
					+ std::to_string(m_timeCreated)
					+ std::to_string(m_timeReceived)
					+ std::to_string(m_receivedFromId)
					+ std::to_string(m_senderId);

	CryptoPP::SHA256 hash;
	byte digest[ CryptoPP::SHA256::DIGESTSIZE ];
	hash.CalculateDigest( digest, reinterpret_cast<byte*>(&msg[0]), msg.length() );

	CryptoPP::HexEncoder encoder;
	std::string output;
	encoder.Attach( new CryptoPP::StringSink( output ) );
	encoder.Put( digest, sizeof(digest) );
	encoder.MessageEnd();

	m_blockHash = output;
	return output;
}

bool Block::VerifyHash(std::string recv_hash) {
	if (m_blockHash.compare(recv_hash) != 0){
		return false; // Different hash received
	}
	return true;
}

void Block::SetTransaction(Transaction t){
	m_transaction = t;
}

int Block::GetBlockHeight(void) const {
	return m_blockHeight;
}

void Block::SetBlockHeight(int height) {
	m_blockHeight = height;
}

int Block::GetContentID(void) const {
	return m_contentID;
}
void Block::SetContentID(int contentID) {
	m_contentID = contentID;
}

int Block::GetDecision(void) const {
	return m_decision;
}
void Block::SetDecision(int decision) {
	m_decision = decision;
}

int Block::GetParentBlockContentID(void) const {
	return m_parentBlockContentID;
}
void Block::SetParentBlockContentID(int parentBlockContentID) {
	m_parentBlockContentID = parentBlockContentID;
}

std::string Block::GetParentBlockHash(void) const{
	return m_parentBlockHash;
}
void Block::SetParentBlockHash(std::string parentBlockHash){
	m_parentBlockHash = parentBlockHash;
}

int Block::GetByteblockSize(void) const {
	return m_byteblockSize;
}

void Block::SetByteblockSize(int byteblockSize){
	m_byteblockSize = byteblockSize;
}

double Block::GetTimeCreated(void) const {
	return m_timeCreated;
}

double Block::GetTimeReceived(void) const {
	return m_timeReceived;
}

Ipv4Address Block::GetReceivedFrom(void) const {
	return m_receivedFrom;
}

void Block::SetReceivedFrom(Ipv4Address receivedFrom) {
	m_receivedFrom = receivedFrom;
}

int Block::GetReceivedFromId(void) const {
	return m_receivedFromId;
}

void Block::SetReceivedFromId(int receivedFrom) {
	m_receivedFromId = receivedFrom;
}

Ipv4Address Block::GetSender(void) const {
	return m_sender;
}

void Block::SetSender(Ipv4Address sender) {
	m_sender = sender;
}

int Block::GetSenderId(void) const {
	return m_senderId;
}

void Block::SetSenderId(int senderId) {
	m_senderId = senderId;
}

bool Block::IsParent(const Block &block) const {
	int blockHeight = GetBlockHeight();
	int contentID = GetContentID();
	if (blockHeight == block.GetBlockHeight() - 1 && contentID == block.GetParentBlockContentID()) {
		return true;
	} else {
		return false;
	}
}


bool Block::IsChild(const Block &block) const {
	int blockHeight = GetBlockHeight();
	int parentContentID = GetParentBlockContentID();
	if (blockHeight == block.GetBlockHeight() + 1 && parentContentID == block.GetContentID()) {
		return true;
	} else {
		return false;
	}
}

void Block::Print(void) {
	std::cout <<"\tTx_hash " << m_transaction.getTxHash() << ",\n" <<
				"\tBlockheight: " << m_blockHeight << ",\n" <<
				"\tContentID: " << m_contentID << ",\n " <<
				"\tParentBlockContentID: " << m_parentBlockContentID << ",\n" <<
				"\tParentBlockHash" << m_parentBlockHash << ",\n" <<
				"\tByteblockSize: " << m_byteblockSize << ",\n" <<
				"\tTimestamp: " << m_timeCreated << ",\n" <<
				"\tReceivedFrom: " << m_receivedFrom << ",\n" <<
				"\tReceivedFromId: " << m_receivedFromId << ",\n" <<
				"\tSender: " << m_sender << ",\n" <<
				"\tSenderId: " << m_senderId << ",\n" <<
				"\tDecision: " << m_decision << 
	std::endl;
}

/*Block& operator= (const Block &blockSrc) {
	m_blockHeight = blockSrc.m_blockHeight;
	m_contentID = blockSrc.m_contentID;
	m_parentBlockContentID = blockSrc.m_parentBlockContentID;
	m_byteblockSize = blockSrc.m_byteblockSize;
	m_timeCreated = blockSrc.m_timeCreated;
	m_timeReceived = blockSrc.m_timeReceived;
	m_receivedFrom = blockSrc.m_receivedFrom;

	return *this;
}*/

bool operator== (const Block &b1, const Block &b2) {
	if (b1.GetBlockHeight() == b2.GetBlockHeight() && b1.GetContentID() == b2.GetContentID()) {
		return true;
	} else {
		return false;
	}
}

std::ostream& operator<< (std::ostream &out, const Block &block) {
	out <<  "(m_blockHeight: " << block.GetBlockHeight() << ", " <<
        "m_contentID: " << block.GetContentID() << ", " <<
        "m_parentBlockContentID: " << block.GetParentBlockContentID() << ", " <<
        "m_parentBlockHash:" << block.GetParentBlockHash() << ", " <<
        "m_byteblockSize: " << block.GetByteblockSize() << ", " <<
        "m_timeCreated: " << block.GetTimeCreated() << ", " <<
        "m_timeReceived: " << block.GetTimeReceived() << ", " <<
        "m_receivedFrom: " << block.GetReceivedFrom() << ", " <<
        "m_receivedFromId: " << block.GetReceivedFromId() <<
        "m_sender: " << block.GetSender() << ", " <<
        "m_senderId: " << block.GetSenderId() <<
        ")";
    return out;
}

class BlockchainChunk : public Block {
public: 
	BlockchainChunk (Transaction transaction, 
		int decision, int blockHeight, int contentId, 
		int chunkID, int parentBlockContentId = 0, 
		std::string parentBlockHash = std::string(""), int blockSizeBytes = 0, 
		double timeCreated = 0, double timeReceived = 0, 
		Ipv4Address receivedFrom = Ipv4Address("0.0.0.0"), int receivedFromId = -1,
		Ipv4Address sender = Ipv4Address("0.0.0.0"), int senderId = -1);

	BlockchainChunk(void);
	BlockchainChunk(const BlockchainChunk &chunkSrc);
	virtual ~BlockchainChunk(void);
	int GetChunkID(void) const;
	void SetChunkID(int chunkID);

	BlockchainChunk& operator= (const BlockchainChunk &chunkSrc);
	friend bool operator== (const BlockchainChunk &chunk, const BlockchainChunk &chunk2);
	friend bool operator< (const BlockchainChunk &chunk, const BlockchainChunk &chunk2);
	friend std::ostream& operator << (std::ostream &out, const BlockchainChunk &chunk);
protected:
	int 			m_chunkID;

};

BlockchainChunk::BlockchainChunk (Transaction transaction, int decision, int blockHeight, 
									int contentID, int chunkID, int parentBlockContentId, 
									std::string parentBlockHash, int blockSizeBytes, 
									double timeCreated, double timeReceived, 
									Ipv4Address receivedFrom, int receivedFromId, 
									Ipv4Address sender, int senderId): 
				Block(transaction, decision, blockHeight, contentID,  
					parentBlockContentId, parentBlockHash, blockSizeBytes, timeCreated, timeReceived, 
					receivedFrom, receivedFromId, sender, senderId) 
{
	m_chunkID = chunkID;
}

BlockchainChunk::BlockchainChunk()
{  
  Transaction t;
  BlockchainChunk(t, -1, -1, -1, -1, -1, std::string(""), -1, -1, -1, Ipv4Address("0.0.0.0"), -1, Ipv4Address("0.0.0.0"), -1);
}

BlockchainChunk::BlockchainChunk (const BlockchainChunk &chunkSrc)
{  
  m_blockHeight = chunkSrc.m_blockHeight;
  m_contentID = chunkSrc.m_contentID;
  m_chunkID = chunkSrc.m_chunkID;
  m_parentBlockContentID = chunkSrc.m_parentBlockContentID;
  m_byteblockSize = chunkSrc.m_byteblockSize;
  m_timeCreated = chunkSrc.m_timeCreated;
  m_timeReceived = chunkSrc.m_timeReceived;
  m_receivedFrom = chunkSrc.m_receivedFrom;
  m_transaction = chunkSrc.m_transaction;
  m_receivedFromId = chunkSrc.m_receivedFromId;
  m_decision = chunkSrc.m_decision;
  m_parentBlockHash = chunkSrc.m_parentBlockHash;
}

BlockchainChunk::~BlockchainChunk (void)
{
}

int BlockchainChunk::GetChunkID (void) const {
	return m_chunkID;
}

void BlockchainChunk::SetChunkID (int chunkID) {
	m_chunkID = chunkID;
}

BlockchainChunk& 
BlockchainChunk::operator= (const BlockchainChunk &chunkSource) {  
	m_blockHeight = chunkSource.m_blockHeight;
	m_contentID = chunkSource.m_contentID;
	m_chunkID = chunkSource.m_chunkID;
	m_parentBlockContentID = chunkSource.m_parentBlockContentID;
	m_byteblockSize = chunkSource.m_byteblockSize;
	m_timeCreated = chunkSource.m_timeCreated;
	m_timeReceived = chunkSource.m_timeReceived;
	m_receivedFrom = chunkSource.m_receivedFrom;
	m_transaction = chunkSource.m_transaction;
	m_receivedFromId = chunkSource.m_receivedFromId;
	m_decision = chunkSource.m_decision;
	m_parentBlockHash = chunkSource.m_parentBlockHash;
	m_senderId = chunkSource.m_senderId;
	m_sender = chunkSource.m_sender;

	return *this;
}

class Blockchain {
public:
	Blockchain (void);

	virtual ~Blockchain (void);

	int GetNoObsoletBlocks (void) const;
	int GetNoOrphan (void) const;
	int GetTotalBlocks (void) const;
	int GetBlockchainHeight (void) const;



	bool HasBlock (const Block &new_block) const;     // Check if a block is in the chain by comparing a received block
	bool HasBlock(int height, int contentID) const;  // Check if a block is in the chain by height and the related content
	bool HasBlock (int height, std::string hash) const;
	bool HasBlock(int contentID, int* decision);

	void ResetChain(void);

	bool HasBlock(int contentID);

	Block ReturnBlock(int height, int contentID);

	bool IsOrphan (const Block &new_block) const;
	bool IsOrphan (int height, int contentID) const;

	bool isObsolet(const Block &new_block) const;

	const Block* GetBlockPointer(const Block& block) const ;

	const std::vector<const Block*> GetChildrenPointerChain(const Block &block);

	const std::vector<const Block*> GetChildrenOrphanPointerChain(const Block &block);


	const Block* GetParent(const Block& block);

	const Block* GetTheTop (void) const;

	void AddBlock(const Block &block);

	void AddOrphan(const Block &block);

	void RemoveOrphan(const Block &block);

	int GetBlocksInFork(void);

	int GetTheLongestFork(void);

	void Print(void);

	friend std::ostream& operator<<(std::ostream &out, const Blockchain &blockchain);
	std::vector<std::vector<Block>> 	m_blocks;
private:
	int 								m_noObsoletBlocks;
	int 								m_totalBlocks;
	
	std::vector<Block> 					m_orphan;
};

Blockchain::Blockchain(void)
{
  m_noObsoletBlocks = 0;
  m_totalBlocks = 0;
  Transaction t;
  Block genesisBlock(t, -1, -1, -1, -1, std::string(""), -1, -1, -1, Ipv4Address("0.0.0.0"), -1, Ipv4Address("0.0.0.0"), -1);
  AddBlock(genesisBlock); 
}

Blockchain::~Blockchain (void)
{
}

int Blockchain::GetNoObsoletBlocks (void) const {
	return m_noObsoletBlocks;
}

int Blockchain::GetNoOrphan (void) const {
	return m_orphan.size();
}
int Blockchain::GetTotalBlocks (void) const {
	return m_totalBlocks;
}

int Blockchain::GetBlockchainHeight (void) const {
	return GetTheTop()->GetBlockHeight();
}

bool Blockchain::HasBlock (const Block &new_block) const {
	if (new_block.GetBlockHeight() > GetTheTop()->GetBlockHeight()) {
		return false;
	} else {
		for (auto const &block: m_blocks[new_block.GetBlockHeight()]) {
			if (block == new_block) {
				return true;
			}
		}
	}
	return false;
}   

bool Blockchain::HasBlock (int height, int contentID) const {
	if (height >  GetTheTop()->GetBlockHeight()){
		return false;
	} else {
		for (auto const  &block:m_blocks[height]){
			if (block.GetBlockHeight() == height && block.GetContentID() == contentID) {
				return true;
			}
		}
	}
	return false;
}

bool Blockchain::HasBlock (int height, std::string hash) const {
	if (height >  GetTheTop()->GetBlockHeight()){
		return false;
	} else {
		for (auto const  &block:m_blocks[height]){
			if (block.GetBlockHash().compare(hash) == 0 ) {
				return true;
			}
		}
	}
	return false;
}

Block Blockchain::ReturnBlock(int height, int contentID) {
	std::vector<Block>::iterator block_it;
	if (height <= GetBlockchainHeight() && height >= 0){
		for (block_it = m_blocks[height].begin(); block_it < m_blocks[height].end(); block_it++) {
			if (block_it->GetBlockHeight() == height && block_it->GetContentID() == contentID) {
				return *block_it;
			}
		}
	}
	for (block_it = m_orphan.begin(); block_it < m_orphan.end(); block_it++){
		if (block_it->GetBlockHeight() == height && block_it->GetContentID() == contentID) {
			return *block_it;
		}
	}
	Transaction t;
	return Block(t, -1, -1, -1, -1, std::string("-1"), -1, -1, -1, Ipv4Address("0.0.0.0"), -1, Ipv4Address("0.0.0.0"), -1);
}

bool Blockchain::IsOrphan (const Block &new_block) const {
	for (auto const &block : m_orphan) {
		if (block == new_block){
			return true;
		}
	}
	return false;
}
bool Blockchain::IsOrphan (int height, int contentID) const {
	for (auto const &block : m_orphan) {
		if (block.GetBlockHeight() == height && block.GetContentID() == contentID){
			return true;
		}
	}
	return false;
}

const Block* Blockchain::GetBlockPointer(const Block &new_block) const {
	for (auto const &block : m_blocks[new_block.GetBlockHeight()]) {
		if (block == new_block) {
			return &block;
		}
	}
	return nullptr;
}

const std::vector<const Block*> Blockchain::GetChildrenPointerChain(const Block &block){
	std::vector <const Block*> children;
	
	int childHeight = block.GetBlockHeight() + 1;

	if (childHeight > GetBlockchainHeight())
		return children;
	for (std::vector<Block>::iterator block_it = m_blocks[childHeight].begin(); block_it != m_blocks[childHeight].end(); block_it++) {
		if(block.IsParent(*block_it)) {
			children.push_back(&(*block_it));
		}
	}
	return children;
}

const std::vector<const Block*> Blockchain::GetChildrenOrphanPointerChain(const Block &block) {
	std::vector<const Block*> children;
	std::vector<Block>::iterator block_it;

	for (block_it = m_orphan.begin(); block_it < m_orphan.end(); block_it++) {
		if (block.IsParent(*block_it)) {
			children.push_back(&(*block_it));
		}
	}
	return children;
}

const Block* Blockchain::GetParent(const Block& block){
	int parentHeight = block.GetBlockHeight() - 1;
	if (parentHeight > GetBlockchainHeight() || parentHeight < 0) {
		return nullptr;
	}
	for (std::vector<Block>::iterator block_it = m_blocks[parentHeight].begin(); block_it != m_blocks[parentHeight].end(); block_it++) {
		if(block.IsChild(*block_it)) {
			return &(*block_it);
		}
	}
	return nullptr;
}

const Block* Blockchain::GetTheTop (void) const {
	return &m_blocks[m_blocks.size() - 1][0];
}

void Blockchain::ResetChain(void) {
	// Para resetar a lockchain, trocar por um vector vazio;
	m_blocks.clear();

	// Fazer as mesmas operações da construção
	m_noObsoletBlocks = 0;
	m_totalBlocks = 0;
	Transaction t;
	Block genesisBlock(t, -1, -1, -1, -1, std::string("-1"), -1, -1, -1, Ipv4Address("0.0.0.0"), -1, Ipv4Address("0.0.0.0"), -1);
	AddBlock(genesisBlock); 
}

void Blockchain::AddBlock(const Block &new_block) {
	if (m_blocks.size() == 0){
		std::vector<Block> newHeight(1, new_block);
		m_blocks.push_back(newHeight);
	} else if (new_block.GetBlockHeight() > GetTheTop()->GetBlockHeight()){
		int rows = new_block.GetBlockHeight() - GetTheTop()->GetBlockHeight() - 1;
		for (int i  = 0; i < rows; i++) {
			std::vector<Block> newHeight;
			m_blocks.push_back(newHeight);
		}
		std::vector<Block> newHeight(1, new_block);
		m_blocks.push_back(newHeight);
	} else {
		if (m_blocks[new_block.GetBlockHeight()].size() > 0) {
			m_noObsoletBlocks++;
		}
		m_blocks[new_block.GetBlockHeight()].push_back(new_block);
	}
	m_totalBlocks++;
}
void Blockchain::AddOrphan(const Block &block) {
	m_orphan.push_back(block);
}

void Blockchain::RemoveOrphan(const Block &block) {
	std::vector<Block>::iterator block_it;
	for(block_it = m_orphan.begin(); block_it != m_orphan.end(); block_it++){		
		if (block == *block_it) {
			break;
		}
	}
	if (block_it == m_orphan.end()) {
		return;
	} else {
		m_orphan.erase(block_it);
	}
}

int Blockchain::GetBlocksInFork(void){
	std::vector<std::vector<Block>>::iterator blockHeight_it;
	int cont = 0;
	for (blockHeight_it = m_blocks.begin(); blockHeight_it < m_blocks.end(); blockHeight_it++) {
		if (blockHeight_it->size() > 1) {
			cont += blockHeight_it->size();
		}
	}
	return cont;
}

int Blockchain::GetTheLongestFork(void) {
	std::vector<std::vector<Block>>::iterator blockHeight_it;
	std::vector<Block>::iterator block_it;
	std::map<int, int> forkedBlocksParentID;
	std::vector<int> newForks;
	int maxSize = 0;

	for (blockHeight_it = m_blocks.begin(); blockHeight_it < m_blocks.end(); blockHeight_it++) {
		if (blockHeight_it->size() > 1 && forkedBlocksParentID.size() == 0) {
			for (block_it = blockHeight_it->begin(); block_it < blockHeight_it->end(); block_it++) {
				forkedBlocksParentID[block_it->GetContentID()] = 1;
			}
		} else if (blockHeight_it->size() > 1) {
			for (block_it = blockHeight_it->begin(); block_it < blockHeight_it->end(); block_it++) {
				std::map<int, int>::iterator mapIndex = forkedBlocksParentID.find(block_it->GetParentBlockContentID());
				if (mapIndex != forkedBlocksParentID.end()) {
					forkedBlocksParentID[block_it->GetContentID()] = mapIndex->second + 1;
					if (block_it->GetContentID() != mapIndex->first) {
						forkedBlocksParentID.erase(mapIndex);
					}
					newForks.push_back(block_it->GetContentID());
				} else {
					forkedBlocksParentID[block_it->GetContentID()] = 1;
				} //end-else
			} //end-for

			for(auto &block : forkedBlocksParentID) {
				if (std::find(newForks.begin(), newForks.end(), block.first) == newForks.end()) {
					if (block.second > maxSize)
						maxSize = block.second;
					forkedBlocksParentID.erase(block.first);
				} // end-if
			} // end-for
		} // end_elseif
		else if (blockHeight_it->size() == 1 && forkedBlocksParentID.size() > 0) {
			for(auto &block : forkedBlocksParentID) {
				if (std::find(newForks.begin(), newForks.end(), block.first) == newForks.end()) {
					if (block.second > maxSize)
						maxSize = block.second;
					forkedBlocksParentID.clear();
				} // end-if
			} // end-for
		} 
	}
	for(auto &block : forkedBlocksParentID) {
		if (std::find(newForks.begin(), newForks.end(), block.first) == newForks.end()) {
			if (block.second > maxSize)
				maxSize = block.second;
		} // end-if
	} // end-for
	return maxSize;
}
void Blockchain::Print(void) {
	int i;
	std::vector< std::vector<Block>>::iterator blockHeight_it;
	std::vector<Block>::iterator  block_it;
	for (blockHeight_it = m_blocks.begin(), i = 0; blockHeight_it < m_blocks.end(); blockHeight_it++, i++) {
		std::cout << "Block height: " << i << std::endl;
		for (block_it = blockHeight_it->begin();  block_it < blockHeight_it->end(); block_it++)
			block_it->Print();
	}
}

bool Blockchain::HasBlock(int contentID) {
	int i;
	std::vector< std::vector<Block>>::iterator blockHeight_it;
	std::vector<Block>::iterator  block_it;
	for (blockHeight_it = m_blocks.begin(), i = 0; blockHeight_it < m_blocks.end(); blockHeight_it++, i++) {
		for (block_it = blockHeight_it->begin();  block_it < blockHeight_it->end(); block_it++){
			if (block_it->GetContentID() == contentID){
				return true;
			}
		}
	}
	return false;
}

bool Blockchain::HasBlock(int contentID, int* decision) {
	int i;
	std::vector< std::vector<Block>>::iterator blockHeight_it;
	std::vector<Block>::iterator  block_it;
	for (blockHeight_it = m_blocks.begin(), i = 0; blockHeight_it < m_blocks.end(); blockHeight_it++, i++) {
		for (block_it = blockHeight_it->begin();  block_it < blockHeight_it->end(); block_it++){
			if (block_it->GetContentID() == contentID){
				*decision = block_it->GetDecision();
				return true;
			}
		}
	}
	return false;
}

bool operator== (const BlockchainChunk &chunk1, const BlockchainChunk &chunk2) {
	if (chunk1.GetBlockHeight() == chunk2.GetBlockHeight() && chunk1.GetContentID() == chunk2.GetContentID() && chunk1.GetChunkID() == chunk2.GetChunkID())
		return true;
	else
		return false;
}

bool operator< (const BlockchainChunk &chunk1, const BlockchainChunk &chunk2)
{
	if (chunk1.GetBlockHeight() < chunk2.GetBlockHeight())
		return true;
	else if (chunk1.GetBlockHeight() == chunk2.GetBlockHeight() && chunk1.GetContentID() < chunk2.GetContentID())
		return true;
	else if (chunk1.GetBlockHeight() == chunk2.GetBlockHeight() && chunk1.GetContentID() == chunk2.GetContentID() && chunk1.GetChunkID() < chunk2.GetChunkID())
		return true;
	else
		return false;
}

std::ostream& operator<< (std::ostream &out, const BlockchainChunk &chunk)
{
	out <<  "(m_blockHeight: " << chunk.GetBlockHeight() << ", " <<
        "m_contentID: " << chunk.GetContentID() << ", " <<
        "m_chunkID: " << chunk.GetChunkID() << ", " <<
        "m_parentBlockContentID: " << chunk.GetParentBlockContentID() << ", " <<
        "m_byteblockSize: " << chunk.GetByteblockSize() << ", " <<
        "m_timeCreated: " << chunk.GetTimeCreated() << ", " <<
        "m_timeReceived: " << chunk.GetTimeReceived() << ", " <<
        "m_receivedFrom: " << chunk.GetReceivedFrom() <<
        "m_receivedFromId: " << chunk.GetReceivedFromId() <<
        "sender: " << chunk.GetReceivedFrom() <<
        "m_receivedFrom: " << chunk.GetReceivedFromId() <<
        ")";

    return out;
}


std::ostream& operator<< (std::ostream &out, Blockchain &blockchain)
{
  
  std::vector< std::vector<Block>>::iterator blockHeight_it;
  std::vector<Block>::iterator  block_it;
  int i;
  
  for (blockHeight_it = blockchain.m_blocks.begin(), i = 0; blockHeight_it < blockchain.m_blocks.end(); blockHeight_it++, i++) 
  {
    out << "  BLOCK HEIGHT " << i << ":\n";
    for (block_it = blockHeight_it->begin();  block_it < blockHeight_it->end(); block_it++)
    {
      out << *block_it << "\n";
    }
  }
  
  return out;
}